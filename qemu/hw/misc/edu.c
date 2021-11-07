/*
 * QEMU qccational PCI device
 *
 * Copyright (c) 2012-2015 Jiri Slaby
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#include "qemu/osdep.h"
#include "qemu/units.h"
#include "hw/pci/pci.h"
#include "hw/hw.h"
#include "hw/pci/msi.h"
#include "qemu/timer.h"
#include "qom/object.h"
#include "qemu/main-loop.h" /* iothread mutex */
#include "qemu/module.h"
#include "qapi/visitor.h"
#include "qc.h"

// complex number operations
Complex conjugate(Complex a) {
	Complex b = {a.real, -1*a.imag};
	return b;
} 

Complex cadd(Complex a, Complex b) {
	Complex sum = {a.real + b.real, a.imag + b.imag};
	return sum;
}

Complex csub(Complex a, Complex b) {
	Complex diff = {a.real - b.real, a.imag - b.imag};
	return diff;
}

Complex cmul(Complex a, Complex b) {
	Complex prod = {a.real * b.real - (a.imag * b.image), a.real * b.imag + (a.imag * b.real)};
	return prod;
}


Complex cdiv(Complex a, Complex b) {
	

}


void vec_push(CVec* vector, Complex el) {
	Item* new_item = (Item*)malloc(sizeof(Item));
//	item->ptr = 

}


#define TYPE_PCI_QC_DEVICE "qc"
typedef struct QCState QCState;
DECLARE_INSTANCE_CHECKER(QCState, QC,
                         TYPE_PCI_QC_DEVICE)

#define FACT_IRQ        0x00000001
#define DMA_IRQ         0x00000100

#define DMA_START       0xFEA40000
#define DMA_SIZE        4096

struct QCState {
    PCIDevice pdev;
    MemoryRegion mmio;

    QemuThread thread;
    QemuMutex thr_mutex;
    QemuCond thr_cond;
    bool stopping;

    uint32_t addr4;
    uint32_t fact;
#define QC_STATUS_COMPUTING    0x01
#define QC_STATUS_IRQFACT      0x80
    uint32_t status;

    uint32_t irq_status;

#define QC_DMA_RUN             0x1
#define QC_DMA_DIR(cmd)        (((cmd) & 0x2) >> 1)
# define QC_DMA_FROM_PCI       0
# define QC_DMA_TO_PCI         1
#define QC_DMA_IRQ             0x4
    struct dma_state {
        dma_addr_t src;
        dma_addr_t dst;
        dma_addr_t cnt;
        dma_addr_t cmd;
    } dma;
    QEMUTimer dma_timer;
    char dma_buf[DMA_SIZE];
    uint64_t dma_mask;
};

static bool qc_msi_enabled(QCState *qc)
{
    return msi_enabled(&qc->pdev);
}

static void qc_raise_irq(QCState *qc, uint32_t val)
{
    qc->irq_status |= val;
    if (qc->irq_status) {
        if (qc_msi_enabled(qc)) {
            msi_notify(&qc->pdev, 0);
        } else {
            pci_set_irq(&qc->pdev, 1);
        }
    }
}

static void qc_lower_irq(QCState *qc, uint32_t val)
{
    qc->irq_status &= ~val;

    if (!qc->irq_status && !qc_msi_enabled(qc)) {
        pci_set_irq(&qc->pdev, 0);
    }
}


static bool within(uint64_t addr, uint64_t start, uint64_t end)
{
    return start <= addr && addr < end;
}

static void qc_check_range(uint64_t addr, uint64_t size1, uint64_t start,
                uint64_t size2)
{
    uint64_t end1 = addr + size1;
    uint64_t end2 = start + size2;

    if (within(addr, start, end2) &&
            end1 > addr && within(end1, start, end2)) {
        return;
    }

    hw_error("QC: DMA range 0x%016"PRIx64"-0x%016"PRIx64
             " out of bounds (0x%016"PRIx64"-0x%016"PRIx64")!",
            addr, end1 - 1, start, end2 - 1);
}

static dma_addr_t qc_clamp_addr(const QCState *qc, dma_addr_t addr)
{
    dma_addr_t res = addr & qc->dma_mask;

    if (addr != res) {
        printf("QC: clamping DMA %#.16"PRIx64" to %#.16"PRIx64"!\n", addr, res);
    }

    return res;
}

static void qc_dma_timer(void *opaque)
{
    QCState *qc = opaque;
    bool raise_irq = false;

    if (!(qc->dma.cmd & QC_DMA_RUN)) {
        return;
    }

    if (QC_DMA_DIR(qc->dma.cmd) == QC_DMA_FROM_PCI) {
       
	uint64_t dst = qc->dma.dst;
	
	qc_check_range(dst, qc->dma.cnt, DMA_START, DMA_SIZE);
        dst -= DMA_START;
	
        pci_dma_read(&qc->pdev, qc_clamp_addr(qc, qc->dma.src),
       		qc->dma_buf + dst, qc->dma.cnt);

    } else {
        uint64_t src = qc->dma.src;
        qc_check_range(src, qc->dma.cnt, DMA_START, DMA_SIZE);
        src -= DMA_START;
        pci_dma_write(&qc->pdev, qc_clamp_addr(qc, qc->dma.dst),
                qc->dma_buf + src, qc->dma.cnt);
    }

    qc->dma.cmd &= ~QC_DMA_RUN;
    if (qc->dma.cmd & QC_DMA_IRQ) {
        raise_irq = true;
    }

    if (raise_irq) {
        qc_raise_irq(qc, DMA_IRQ);
    }
}



static void dma_rw(QCState *qc, bool write, dma_addr_t *val, dma_addr_t *dma,
                bool timer)
{
    if (write && (qc->dma.cmd & QC_DMA_RUN)) {
        return;
    }

    if (write) {
        *dma = *val;
    } else {
        *val = *dma;
    }

    if (timer) {
        timer_mod(&qc->dma_timer, qemu_clock_get_ms(QEMU_CLOCK_VIRTUAL) + 100);
    }
}



static uint64_t qc_mmio_read(void *opaque, hwaddr addr, unsigned size)
{
    QCState *qc = opaque;
    uint64_t val = ~0ULL;

    if (addr < 0x80 && size != 4) {
        return val;
    }

    if (addr >= 0x80 && size != 4 && size != 8) {
        return val;
    }

    switch (addr) {
    case 0x00:
        val = 0x5143534D;
	
        break;
    case 0x04:
        val = qc->addr4;
        break;
    case 0x08:
        qemu_mutex_lock(&qc->thr_mutex);
        val = qc->fact;
        qemu_mutex_unlock(&qc->thr_mutex);
        break;
    case 0x20:
        val = qatomic_read(&qc->status);
        break;
    case 0x24:
        val = qc->irq_status;
        break;
    case 0x80:
        dma_rw(qc, false, &val, &qc->dma.src, false);
        break;
    case 0x88:
        dma_rw(qc, false, &val, &qc->dma.dst, false);
        break;
    case 0x90:
        dma_rw(qc, false, &val, &qc->dma.cnt, false);
        break;
    case 0x98:
        dma_rw(qc, false, &val, &qc->dma.cmd, false);
        break;
    }

    return val;
}

static void qc_mmio_write(void *opaque, hwaddr addr, uint64_t val,
                unsigned size)
{
    QCState *qc = opaque;

    if (addr < 0x80 && size != 4) {
        return;
    }

    if (addr >= 0x80 && size != 4 && size != 8) {
        return;
    }

    switch (addr) {
    case 0x04:
        qc->addr4 = ~val;
        break;
    case 0x08:
        if (qatomic_read(&qc->status) & QC_STATUS_COMPUTING) {
            break;
        }
        /* QC_STATUS_COMPUTING cannot go 0->1 concurrently, because it is only
         * set in this function and it is under the iothread mutex.
         */
        qemu_mutex_lock(&qc->thr_mutex);
        qc->fact = val;
        qatomic_or(&qc->status, QC_STATUS_COMPUTING);
        qemu_cond_signal(&qc->thr_cond);
        qemu_mutex_unlock(&qc->thr_mutex);
        break;
    case 0x20:
        if (val & QC_STATUS_IRQFACT) {
            qatomic_or(&qc->status, QC_STATUS_IRQFACT);
        } else {
            qatomic_and(&qc->status, ~QC_STATUS_IRQFACT);
        }
        break;
    case 0x60:
        qc_raise_irq(qc, val);
        break;
    case 0x64:
        qc_lower_irq(qc, val);
        break;
    case 0x80:
        dma_rw(qc, true, &val, &qc->dma.src, false);
        break;
    case 0x88:
        dma_rw(qc, true, &val, &qc->dma.dst, false);
        break;
    case 0x90:
        dma_rw(qc, true, &val, &qc->dma.cnt, false);
        break;
    case 0x98:
        if (!(val & QC_DMA_RUN)) {
            break;
        }
        dma_rw(qc, true, &val, &qc->dma.cmd, true);
        break;
    }
}


static const MemoryRegionOps qc_mmio_ops = {
    .read = qc_mmio_read,
    .write = qc_mmio_write,
    .endianness = DEVICE_NATIVE_ENDIAN,
    .valid = {
        .min_access_size = 4,
        .max_access_size = 8,
    },
    .impl = {
        .min_access_size = 4,
        .max_access_size = 8,
    },

};

/*
 * We purposely use a thread, so that users are forced to wait for the status
 * register.
 */
static void *qc_fact_thread(void *opaque)
{
    QCState *qc = opaque;

    while (1) {
        uint32_t val, ret = 1;

        qemu_mutex_lock(&qc->thr_mutex);
        while ((qatomic_read(&qc->status) & QC_STATUS_COMPUTING) == 0 &&
                        !qc->stopping) {
            qemu_cond_wait(&qc->thr_cond, &qc->thr_mutex);
        }

        if (qc->stopping) {
            qemu_mutex_unlock(&qc->thr_mutex);
            break;
        }

        val = qc->fact;
        qemu_mutex_unlock(&qc->thr_mutex);

        while (val > 0) {
            ret *= val--;
        }

        /*
         * We should sleep for a random period here, so that students are
         * forced to check the status properly.
         */

        qemu_mutex_lock(&qc->thr_mutex);
        qc->fact = ret;
        qemu_mutex_unlock(&qc->thr_mutex);
        qatomic_and(&qc->status, ~QC_STATUS_COMPUTING);

        if (qatomic_read(&qc->status) & QC_STATUS_IRQFACT) {
            qemu_mutex_lock_iothread();
            qc_raise_irq(qc, FACT_IRQ);
            qemu_mutex_unlock_iothread();
        }
    }

    return NULL;
}

static void pci_qc_realize(PCIDevice *pdev, Error **errp)
{
    QCState *qc = QC(pdev);

    uint8_t *pci_conf = pdev->config;

    pci_config_set_interrupt_pin(pci_conf, 1);

    if (msi_init(pdev, 0, 1, true, false, errp)) {
        return;
    }

    timer_init_ms(&qc->dma_timer, QEMU_CLOCK_VIRTUAL, qc_dma_timer, qc);

    qemu_mutex_init(&qc->thr_mutex);
    qemu_cond_init(&qc->thr_cond);
    qemu_thread_create(&qc->thread, "qc", qc_fact_thread,
                       qc, QEMU_THREAD_JOINABLE);

    memory_region_init_io(&qc->mmio, OBJECT(qc), &qc_mmio_ops, qc,
                    "qc-mmio", 1 * MiB);
    pci_register_bar(pdev, 0, PCI_BASE_ADDRESS_SPACE_MEMORY, &qc->mmio);
}

static void pci_qc_uninit(PCIDevice *pdev)
{
    QCState *qc = QC(pdev);

    qemu_mutex_lock(&qc->thr_mutex);
    qc->stopping = true;
    qemu_mutex_unlock(&qc->thr_mutex);
    qemu_cond_signal(&qc->thr_cond);
    qemu_thread_join(&qc->thread);

    qemu_cond_destroy(&qc->thr_cond);
    qemu_mutex_destroy(&qc->thr_mutex);

    timer_del(&qc->dma_timer);
    msi_uninit(pdev);
}

static void qc_instance_init(Object *obj)
{
    QCState *qc = QC(obj);

    qc->dma_mask = (1UL << 28) - 1;
    object_property_add_uint64_ptr(obj, "dma_mask",
                                   &qc->dma_mask, OBJ_PROP_FLAG_READWRITE);
}

static void qc_class_init(ObjectClass *class, void *data)
{
    DeviceClass *dc = DEVICE_CLASS(class);
    PCIDeviceClass *k = PCI_DEVICE_CLASS(class);

    k->realize = pci_qc_realize;
    k->exit = pci_qc_uninit;
    k->vendor_id = 0xC0DE/*PCI_VENDOR_ID_QEMU*/;
    k->device_id = 0xDEAD/*0x11e8*/; 
    k->revision = 0x10;
    k->class_id = PCI_CLASS_OTHERS;
    set_bit(DEVICE_CATEGORY_MISC, dc->categories);

}

static void pci_qc_register_types(void)
{
    static InterfaceInfo interfaces[] = {
        { INTERFACE_CONVENTIONAL_PCI_DEVICE },
        { },
    };
    static const TypeInfo qc_info = {
        .name          = TYPE_PCI_QC_DEVICE,
        .parent        = TYPE_PCI_DEVICE,
        .instance_size = sizeof(QCState),
        .instance_init = qc_instance_init,
        .class_init    = qc_class_init,
        .interfaces = interfaces,
    };

    type_register_static(&qc_info);
}

type_init(pci_qc_register_types)
