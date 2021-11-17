#ifndef QC_H
#define QC_H

#include <stddef.h>
#include <math.h>
#include "qc.h"

#define CVAL(ADDR) (*(Complex*)(ADDR))
#define CPTR(ADDR) ((Complex*)(ADDR))
#define VVAL(ADDR) (*(CVec*)(ADDR))
#define VPTR(ADDR) ((CVec*)(ADDR))


typedef struct Complex {
	double real;
	double imag;

} Complex;


Complex cadd(Complex a, Complex b);

Complex csub(Complex a, Complex b);

Complex cmul(Complex a, Complex b);

Complex cdiv(Complex a, Complex b);


typedef struct Item {
	void* val;
	struct Item* prev;
	struct Item* next;

} Item;


typedef struct CVec {
	unsigned int size;
	Item* head;
	Item* tail;
		

} CVec;


typedef struct Matrix {
	unsigned int m;
	unsigned int n;
	Item* head;
	Item* tail;

} Matrix;



void vec_push(CVec* vector, void* el);

CVec* create_vec();



// get element at specified index
Item* subscript(CVec* v, unsigned int index);

// get complex number at index
Complex cindex(CVec* v, unsigned int index);


// get vector at index (presumably in a matrix)
CVec vindex(CVec* m, unsigned int index);



// inner product (dot for short)
Complex dot(CVec* a, CVec* b);


// check if vectors are orthogonal
int is_orthogonal(CVec* a, CVec* b);


// print a formatted complex number
void printc(Complex a);


// apply a linear operation (usually matrix in a basis) to a vector
CVec* linop(CVec* op, CVec* state);


// absolute square
Complex absq(Complex a);


// print probabilities
void print_prob(CVec* v);


// print statevector
void print_vector(CVec* v);

// get the transpose of a matrix (switch rows + columns)
CVec* transpose(CVec * matrix);


 
// calculate tensor product for two vectors
CVec* tensor_vec(CVec* a, CVec* b);



// calculate tensor for matrices
CVec* tensor(CVec* a, CVec* b);


void print_matrix(CVec* m);

Complex* create_complex(double real, double imag);

CVec* create_matrix(Complex* a, Complex* b, Complex* c, Complex* d);


CVec* cnot_gate(void);


CVec* h_gate(void);


CVec* id_gate(void);


CVec* x_gate(void);


CVec* z_gate(void);


CVec* y_gate(void);



#endif
