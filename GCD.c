#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

unsigned int GCD(unsigned int a, unsigned int b)
{
	unsigned int temp=0;
	if(a==0)
		return b;
	if(b==0)
		return a;
	while(b!=0)
	{
		temp=b;
		b=a%b;
		a=temp;
	}
	return a;
}

int main()
{
	unsigned int a=0, b=0, gcd=0;
	double time=0.0;
	printf("Enter first number: ");
	scanf("%u", &a);
	printf("Enter second number: ");
	scanf("%u", &b);
	clock_t start = clock();
	gcd = GCD(a,b);
	clock_t end = clock();
	time= ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("GCD = %u\n",gcd);
	printf("Time taken: %f seconds\n", time);
	return 0;
}

//calculating time has been referred from GeeksForGeeks
