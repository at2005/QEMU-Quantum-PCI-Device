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



void linop(CVec* op, CVec* state) {
	
	for(int i = 0; i < op->size; i++) {
		CVec vec = vindex(op,i);
		printc(dot(&vec, state));
	
	}


}


int main() {
	Complex test1 = {1/sqrt(2),0};
	Complex test2 = {1/sqrt(2),0};
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
	
	linop(m, v);
		

	/*vec_push(v, &test1);
	vec_push(v, &test2);
	vec_push(v, &test6);
		
	Complex test3 = {1,0};	
	Complex test4 = {4,0};	
	//Complex test5 = {7,0};	

	CVec* v2 = create_vec();	
	vec_push(v2, &test3);	
	vec_push(v2, &test4);	
	vec_push(v2, &test5);
	
	printc(dot(v,v2));	

	//printc(CVAL(subscript(v,0)->val));
	*/
	return 0;


}
