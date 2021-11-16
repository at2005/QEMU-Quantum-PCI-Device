#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "qc.h"

#define CVAL(ADDR) (*(Complex*)(ADDR))
#define CPTR(ADDR) ((Complex*)(ADDR))
#define VVAL(ADDR) (*(CVec*)(ADDR))
#define VPTR(ADDR) ((CVec*)(ADDR))

// complex number operations
Complex conjugate(Complex a) {
	Complex b = {a.real, -1*a.imag};
	return b;
} 


// addition
Complex cadd(Complex a, Complex b) {
	Complex sum = {a.real + b.real, a.imag + b.imag};
	return sum;
}


// subtraction
Complex csub(Complex a, Complex b) {
	Complex diff = {a.real - b.real, a.imag - b.imag};
	return diff;
}

// multiplication 
Complex cmul(Complex a, Complex b) {
	Complex prod = {a.real * b.real - (a.imag * b.imag), a.real * b.imag + (a.imag * b.real)};
	return prod;
}


// division (multiplication by conjugate transpose)
Complex cdiv(Complex a, Complex b) {
	Complex ratio;
	ratio.real = a.real * b.real + b.imag * a.imag;
	ratio.real /= b.real * b.real + b.imag * b.imag;
	ratio.imag = (b.real * a.imag) - a.real * b.imag;
	ratio.imag /= b.real * b.real + b.imag * b.imag;
	return ratio;
}


// add element to vector
void vec_push(CVec* vector, void* el) {
	// allocate and assign new item
	Item* new_item = (Item*)malloc(sizeof(Item));
	new_item->val = el;
	
	
	// check if vector size is zero so that we don't dereference any nullptr
	if(vector->size != 0) {
		vector->tail->next = new_item;
		new_item->prev = vector->tail;
		vector->tail = new_item;
		new_item->next = vector->head;
		vector->head->prev = new_item;
		vector->size++;
		return;
	}
	
	// case for first element, set up vector head and tail
	vector->head = new_item;
	vector->tail = new_item;
	new_item->prev = new_item;
	new_item->next = new_item;
	vector->size++;
	
}


// create a vector
CVec* create_vec() {
	// allocate vector and initialize
	CVec* new_vec = (CVec*)malloc(sizeof(CVec));
	new_vec->head = (Item*)NULL;
	new_vec->tail = (Item*)NULL;
	new_vec->size = 0;
	return new_vec;
	

}


// get element at specified index
Item* subscript(CVec* v, unsigned int index) {
	Item* init = v->head;
	for(int i = 0; i < index; i++) {
		init = init->next;
	
	}

	return init;

}

// get complex number at index
Complex cindex(CVec* v, unsigned int index) {
	Item* it = subscript(v, index);
	return *(Complex*)(it->val);
	
}


// get vector at index (presumably in a matrix)
CVec vindex(CVec* m, unsigned int index) {
	Item* it = subscript(m,index);
	return *(CVec*)(it->val);


}


// inner product (dot for short)
Complex dot(CVec* a, CVec* b) {
	Complex total = {0,0};
	for(int i = 0; i < a->size; i++) {
		// we take the conjugate of the bra because that is the definition of inner product
		total = cadd(total, cmul(conjugate(CVAL(subscript(a,i)->val)), CVAL(subscript(b,i)->val)));
	
	}	
	
	return total;
}


// check if vectors are orthogonal
int is_orthogonal(CVec* a, CVec* b) {
	if(dot(a,b).real == 0 && dot(a,b).imag == 0) {
		return 1;
	
	}

	return 0;

}


// print a formatted complex number
void printc(Complex a) {
	printf("%f+%fi\n", a.real, a.imag);

}


// apply a linear operation (usually matrix in a basis) to a vector
CVec* linop(CVec* op, CVec* state) {
	CVec* res = create_vec();
	for(int i = 0; i < op->size; i++) {
		// get row at specified index
		CVec vec = vindex(op,i);
		Complex* mul_res = (Complex*)malloc(sizeof(Complex));
		// get the inner product of row + column vector
		*mul_res = dot(&vec, state);
		// append to result vector
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


// print probabilities
void print_prob(CVec* v) {
	for(int i = 0; i < v->size; i++) {
		printf("P(%d) is %f% \n", i, absq(CVAL(subscript(v,i)->val)).real * 100);
	
	}


}


// print statevector
void print_vector(CVec* v) {
	for(int i = 0; i < v->size; i++) {
		Complex num = CVAL(subscript(v,i)->val);
		printf("%f+%fi|%d>\n",num.real, num.imag,i);
	//	if(i+1 < v->size) printf("+ ");
	
	}

//	printf("\n");



}

// get the transpose of a matrix (switch rows + columns)
CVec* transpose(CVec * matrix) {
	CVec* new_matrix = create_vec();
	// iterate multiple times because we need to extract elements from the column vectors but 
	// matrix is defined as a vector of row-vectors, so we have to iterate b^2 times where b is 
	// the size of the matrix
	for(int i = 0; i < matrix->size; i++) {
		CVec* col = create_vec();
		for(int j = 0; j < matrix->size; j++) {
			// get the vector at index j and store in v
			CVec* v = VPTR(subscript(matrix,j)->val);
			// get the element at index i in v
			Complex* c = CPTR(subscript(v,i)->val);
			// append vector to matrix
			vec_push(col,c);
		
		}

		vec_push(new_matrix,col);	
	
	}

	return new_matrix;

} 


 
// calculate tensor product for two vectors
CVec* tensor_vec(CVec* a, CVec* b) {
	CVec* new_vec = create_vec();
	// iterate over each vector and each vector's elements and calculate tensor product
	// e.g. (a_0*b_0, a_0*b_1, a_1*b0, a1*b1 ... a_n*b_n)
	for(int i = 0; i < a->size; i++) {
		Complex c1 = CVAL(subscript(a,i)->val);
		for(int j = 0; j < b->size; j++) {
			Complex c2 = CVAL(subscript(b,j)->val);
			Complex* res = (Complex*)malloc(sizeof(Complex));
			*res = cmul(c1,c2);
			vec_push(new_vec, res);

		}
	
	}	

	return new_vec;
}



// calculate tensor for matrices
CVec* tensor(CVec* a, CVec* b) {
	CVec* tensor_prod = create_vec();
	// transpose matrix b so that we can easily access columns as rows
	CVec* b_dg = transpose(b);
	// iterate and calculate vector tensor product and add to matrix
	for(int i = 0; i < a->size; i++) {
		CVec* col = VPTR(subscript(a,i)->val);
		for(int j = 0; j < b_dg->size; j++) {
			CVec* row = VPTR(subscript(b,j)->val);
			vec_push(tensor_prod, tensor_vec(col,row)); 	

		}
	
	}

	for(int i = 0; i < b_dg->size; i++) {
		free(VPTR(subscript(b_dg,i)->val));
	}

	free(b_dg);
	
	return tensor_prod;	
}


void print_matrix(CVec* m) {
	for(int i = 0; i < m->size; i++) {
		CVec* v = VPTR(subscript(m,i)->val);
		for(int j = 0; j < v->size; j++) {
			Complex c = CVAL(subscript(v,j)->val);
			printf("%f + %fi  ", c.real, c.imag);
			
		}

		printf("\n");
	}

}

Complex* create_complex(double real, double imag) {
	Complex* a = (Complex*)malloc(sizeof(Complex));
	a->real = real;
	a->imag = imag;
	return a;

}

CVec* create_matrix(Complex* a, Complex* b, Complex* c, Complex* d) {
	CVec* matrix = create_vec();	
	CVec* row1 = create_vec();
	CVec* row2 = create_vec();

	vec_push(row1, a);
	vec_push(row1, b);
	vec_push(row2, c);
	vec_push(row2, d);

	vec_push(matrix, row1);
	vec_push(matrix, row2);
	
	return matrix;

}


CVec* cnot_gate() {
	CVec* matrix = create_vec();
	CVec* row1 = create_vec();
	CVec* row2 = create_vec();
	CVec* row3 = create_vec();
	CVec* row4 = create_vec();
	
	Complex* m1 = create_complex(1,0);
	Complex* m2 = create_complex(0,0);
	Complex* m3 = create_complex(0,0);
	Complex* m4 = create_complex(0,0);
	
	vec_push(row1, m1);
	vec_push(row1, m2);
	vec_push(row1, m3);
	vec_push(row1, m4);

	Complex* m5 = create_complex(0,0);
	Complex* m6 = create_complex(1,0);
	Complex* m7 = create_complex(0,0);
	Complex* m8 = create_complex(0,0);

	vec_push(row2, m5);
	vec_push(row2, m6);
	vec_push(row2, m7);
	vec_push(row2, m8);

	Complex* m9 = create_complex(0,0);
	Complex* m10 = create_complex(0,0);
	Complex* m11 = create_complex(0,0);
	Complex* m12 = create_complex(1,0);
	
	vec_push(row3, m9);
	vec_push(row3, m10);
	vec_push(row3, m11);
	vec_push(row3, m12);

	Complex* m13 = create_complex(0,0);
	Complex* m14 = create_complex(0,0);
	Complex* m15 = create_complex(1,0);
	Complex* m16 = create_complex(0,0);

	vec_push(row4, m13);
	vec_push(row4, m14);
	vec_push(row4, m15);
	vec_push(row4, m16);

	vec_push(matrix, row1);
	vec_push(matrix, row2);
	vec_push(matrix, row3);
	vec_push(matrix, row4);
	
	return matrix;
}


CVec* h_gate() {

	Complex* m1 = create_complex(1/sqrt(2), 0);
	Complex* m2 = create_complex(1/sqrt(2), 0);
	Complex* m3 = create_complex(1/sqrt(2), 0);
	Complex* m4 = create_complex(-1/sqrt(2), 0);


	CVec* hadamard_2d = create_matrix(m1,m2,m3,m4);

	return hadamard_2d;
	

}


CVec* id_gate() {
	Complex* m1 = create_complex(1,0);
	Complex* m2 = create_complex(0,0);
	Complex* m3 = create_complex(0,0);
	Complex* m4 = create_complex(1,0);

	CVec* id_2d = create_matrix(m1, m2, m3, m4);
	return id_2d;

}


CVec* x_gate() {
	Complex* m1 = create_complex(0,0);
	Complex* m2 = create_complex(1,0);
	Complex* m3 = create_complex(1,0);
	Complex* m4 = create_complex(0,0);

	CVec* x_2d = create_matrix(m1, m2, m3, m4);
	return x_2d;

}


CVec* z_gate() {
	Complex* m1 = create_complex(1,0);
	Complex* m2 = create_complex(0,0);
	Complex* m3 = create_complex(0,0);
	Complex* m4 = create_complex(-1,0);

	CVec* z_2d = create_matrix(m1, m2, m3, m4);
	return z_2d;


}


CVec* y_gate() {

	Complex* m1 = create_complex(0,0);
	Complex* m2 = create_complex(0,-1);
	Complex* m3 = create_complex(0,1);
	Complex* m4 = create_complex(0,0);

	CVec* y_2d = create_matrix(m1, m2, m3, m4);
	return y_2d;


}



int main() {
	Complex test1 = {1/sqrt(2),0};
	Complex test2 = {1/sqrt(2),0};
	//Complex test6 = {12,0};
	CVec* v = create_vec();
	vec_push(v,&test1);
	vec_push(v,&test2);

	Complex test3 = {1,0};
	Complex test4 = {0,0};
	CVec* v2 = create_vec();
	vec_push(v2, &test3);
	vec_push(v2, &test4);


	CVec* d4 = tensor_vec(v,v2);
	

	//CVec* h = h_gate();
	//CVec* h4 = tensor(h,h);
	
	CVec* cnot = cnot_gate();

	d4 = linop(cnot, d4);
	print_prob(d4);
		

	return 0;


}
