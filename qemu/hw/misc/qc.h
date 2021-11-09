
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
