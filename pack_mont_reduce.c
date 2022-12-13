
#include "finitefield.h"
#include "extensions.h"

// this package is all about the reduction is base on the montgomery reduction
extern bigInteger p;
extern bigInteger b;



void bigIntVagueReduction(bigInteger *a, const bigInteger *p){

   if ( bigIntSignCom(a, p) == LESS &&  a->sign == POSITIVE){

      return;
      
   }

   if (a->sign == POSITIVE){
   
   int divbit = bigIntCountBits(a);
   //int prime = bigIntCountBits(p);
   int k = divbit - 448;   // p is a 447bit number
    
   bigInteger new;
   bigIntInitial(&new);
   //int count = 0;

   while( k >=3 ){

      bigIntMul2k(p,k,&new);

      bigIntSubtraction(a,&new,a);

      //count++;
      //printf("the count number = %d \n", count);

      divbit = bigIntCountBits(a);

      k = divbit - 448;

   }

   while( (bigIntSignCom(a, p) == EQUAL)||(bigIntSignCom(a, p) == GREATER)){
      bigIntSubtraction(a,p,a);  //a is the final value
   }

   return;

}

if (a->sign == NEGATIVE){

   int divbit = bigIntCountBits(a);
   //int prime = bigIntCountBits(p);
   int k = divbit - 448;   // p is a 447bit number

   bigInteger new;
   bigIntInitial(&new);
   //int count = 0;

   while( k >=3 ){

      bigIntMul2k(p,k,&new);

      bigIntAddition(a,&new,a);

      //count++;
      //printf("the count number = %d \n", count);

      divbit = bigIntCountBits(a);

      k = divbit - 448;

   }

   while( (bigIntUnCom(a,p) == EQUAL)||(bigIntUnCom(a, p) == GREATER) || a->sign == NEGATIVE){
      bigIntAddition(a,p,a);  //a is the final value
   }

   return;

}

}




/* Montgomery refers to libtommath*/

void bigIntGenMontReduce(bigInteger *x){

	bigIntDigit su = 4681u;
   passModPrime();
	bigIntMontReduce(x,&p,su);

	generMul(x,&b,x);
  // printf("After Mont, we multi with the b, THE RESULT IS !!! \n");
  // draw(x);

   //bigIntVagueReduction(x,&p);
	bigIntGenMod(x,&p,x);

}

/* su is the setup parameter, which is precomputed according to the modulus p
   su = 4681u;(unsigned int_32)
   and x is the product of two elements in Fp,hence, if x > p, p < x <= (p-1)^2
*/
bigIntCatch bigIntMontReduce(bigInteger *x, const bigInteger *p, bigIntDigit su)
{
   bigIntCatch err; // catch the error
   int xi, newoccu;

   newoccu = (p->occupied * 2) + 6;
   
   /* increase the digit of the input, we will add p's multiple on it*/
   if ((err = bigIntIncrease(x, newoccu)) != runSuccess) {
      return err;
   }
   x->occupied = newoccu;

   for (xi = 0; xi < p->occupied; xi++) {
      int yi; // yi use for loop, for compute every iteration of the result of new x  
      /*
      we have computed the  value of -1/p mod R, is su
      every loop we should pre-compute the value of miumiu =  -xi/p mod R = xi*su mod R
      xi is [0, 2^15-1], su = 4681u, it does not occur overflow!
      */
      bigIntDigit miumiu = ((x->adp[xi] * su) & ANDR_1); 

      bigIntDigit carry = 0;
      /* this "for yi loop" is to compute the every "loop of xi": the result of x, 
      x = x + miumiu * p * R^i */
      for (yi = 0; yi < p->occupied; yi++) {
         /* 
         because in every xi loop, the value of x would continue update
         we donot need to care about the digit of x less than xi
         hence the accumulation is base on x->adp[xi + yi], yi is start from 0
          */
         bigIntDigit ua = (miumiu * p->adp[yi] + carry + x->adp[xi + yi]);// sum of

         /* carry, pass to the next digit*/
         carry = (bigIntDigit)(ua >> DefaultBitOfDigit);

         /*fix the digit of result of x of in every "for xi loop" */
         x->adp[xi + yi] = (ua & ANDR_1);// AND R-1 operation is equal to mod R
      }
      /* pass up carries, if carry is not 0*/
      while (carry != 0u) {
         x->adp[xi + yi] += carry;
         carry = x->adp[xi + yi] >> DefaultBitOfDigit; // pass the carry to next digit
         x->adp[xi + yi++] &= ANDR_1; // fix the digit
      }
   }

   /* now, the digit from 0 to p.occupied of x is 0, we need to right shift p.occupied
      the result of x mod p is equivalent to x/R**(-p.occupied) mod p
   */
   bigIntRightShift(x, p->occupied);
   bigIntAdjust(x);

   return runSuccess;
}








// discard the following code












// the below function is used to check our setup value is 4681u or not

// bigIntDigit setup = 4681u;


/* setups the montgomery reduction stuff */
bigIntCatch mp_montgomery_setup(const bigInteger *n, bigIntDigit *rho)
{
   bigIntDigit x, b;


   /* fast inversion mod 2**k
    *
    * Based on the fact that
    *
    * XA = 1 (mod 2**n)  =>  (X(2-XA)) A = 1 (mod 2**2n)
    *                    =>  2*X*A - X*X*A*A = 1
    *                    =>  2*(1) - (1)     = 1
    */
   b = n->adp[0];

   if ((b & 1u) == 0u) {
      return runInvalid;
   }

   x = (((b + 2u) & 4u) << 1) + b; /* here x*a==1 mod 2**4 */
   x *= 2u - (b * x);              /* here x*a==1 mod 2**8 */
   x *= 2u - (b * x);              /* here x*a==1 mod 2**16 */


   /* rho = -1/m mod b */
   *rho = (bigIntDigit)(((uSignInt)1 << (uSignInt)DefaultBitOfDigit) - x) & ANDR_1;
   

   return runSuccess;
}


