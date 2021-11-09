#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "qc.h"

#define CVAL(ADDR) (*(Complex*)(ADDR))
#define VVAL(ADDR) (*(CVec*)(ADDR))
#define VPTR(ADDR) ((CVec*)(ADDR))

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
	Complex prod = {a.real * b.real - (a.imag * b.imag), a.real * b.imag + (a.imag * b.real)};
	return prod;
}


Complex cdiv(Complex a, Complex b) {
	Complex ratio;
	ratio.real = a.real * b.real + b.imag * a.imag;
	ratio.real /= b.real * b.real + b.imag * b.imag;
	ratio.imag = (b.real * a.imag) - a.real * b.imag;
	ratio.imag /= b.real * b.real + b.imag * b.imag;
	return ratio;
}


void vec_push(CVec* vector, void* el) {
	Item* new_item = (Item*)malloc(sizeof(Item));
	new_item->val = el;
	
	
	if(vector->size != 0) {
		vector->tail->next = new_item;
		new_item->prev = vector->tail;
		vector->tail = new_item;
		new_item->next = vector->head;
		vector->head->prev = new_item;
		vector->size++;
		return;
	}
	

	vector->head = new_item;
	vector->tail = new_item;
	new_item->prev = new_item;
	new_item->next = new_item;
	vector->size++;
	
}



CVec* create_vec() {
	CVec* new_vec = (CVec*)malloc(sizeof(CVec));
	new_vec->head = (Item*)NULL;
	new_vec->tail = (Item*)NULL;
	new_vec->size = 0;
	return new_vec;
	

}


Item* subscript(CVec* v, unsigned int index) {
	Item* init = v->head;
	for(int i = 0; i < index; i++) {
		init = init->next;
	
	}

	return init;

}

Complex cindex(CVec* v, unsigned int index) {
	Item* it = subscript(v, index);
	return *(Complex*)(it->val);
	
}


Complex dot(CVec* a, CVec* b) {
	Complex total = {0,0};
	for(int i = 0; i < a->size; i++) {
		total = cadd(total, cmul(CVAL(subscript(a,i)->val), CVAL(subscript(b,i)->val)));
	
	}	
	
	return total;
}


int is_orthogonal(CVec* a, CVec* b) {
	if(dot(a,b).real == 0 && dot(a,b).imag == 0) {
		return 1;
	
	}

	return 0;

}



void printc(Complex a) {
	printf("%f+%fi\n", a.real, a.imag);

}


CVec vindex(CVec* m, unsigned int index) {
	Item* it = subscript(m,index);
	return *(CVec*)(it->val);


}



CVec* linop(CVec* op, CVec* state) {
	CVec* res = create_vec();
	for(int i = 0; i < op->size; i++) {
		CVec vec = vindex(op,i);
		Complex* mul_res = (Complex*)malloc(sizeof(Complex));
		*mul_res = dot(&vec, state);
		vec_push(res, mul_res);
	
	}
	
	return res;


}


// absolute square
Complex absq(Complex a) {
	Complex b = a;
	b.imag *= -1;
	return cmul(a,b);


}


void print_prob(CVec* v) {
	for(int i = 0; i < v->size; i++) {
		printf("P(%d) is %f% \n", i, absq(CVAL(subscript(v,i)->val)).real * 100);
	
	}


}

int main() {
	Complex test1 = {sin(2),0};
	Complex test2 = {cos(2),0};
	//Complex test6 = {12,0};
	
	CVec* v = create_vec();
	vec_push(v,&test1);
	vec_push(v,&test2);

	
	Complex m1 = {1/sqrt(2),0};
	Complex m2 = {1/sqrt(2),0};
	Complex m3 = {1/sqrt(2),0};
	Complex m4 = {-1/sqrt(2),0};
	
	CVec* c1 = create_vec();
	
	vec_push(c1, &m1);
	vec_push(c1, &m2);
	
	CVec* c2 = create_vec();

	vec_push(c2, &m3);
	vec_push(c2, &m4);
	
	CVec* m = create_vec();
	vec_push(m, c1);
	vec_push(m, c2);
	
	v = linop(m,v);
	print_prob(v);

	v = linop(m,v);
	print_prob(v);
	

	return 0;


}
