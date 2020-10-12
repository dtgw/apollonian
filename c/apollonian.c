/* ***********************************************************
 * Implementation of Apollonian Cell Encoders
 * Initial demonstration version.
 * Authors: Thomas Given-Wilson
 * Version: 0.1.1
 * Licence: 
 * Notes: https://hal.inria.fr/hal-01571226
 * */


// Included for standard libraries, required
#include <stdlib.h>

// Included for numerical types, see notes on number type below
#include <stdint.h>

// Included for CHAR_BIT to be absolutely sure
#include <limits.h>

// Included for seeding random, can be replaced if random is changed
#include <time.h>

// Included for printing, not required if printing (debugging) removed
#include <stdio.h>


// Below for prettiness
#define TRUE 1
#define FALSE 0

// Define the number type, here uint64_t for illustration
// Note that several operators must be defined for this type
// The list of operators is below (also many of these could be changed/reduced):
// Comparison: == != < <= > >= 
// Arithmetic: + - * / % ++ <<
// Other: =
#define NUM_TYPE uint64_t
#define NUM_FMT  "%lu"
// TODO: change the above to something LARGER


/* Debug level
 * Can be used for different levels of debugging.
 * Proposed (not necessarily implemented level below):
 * 0 - no output
 * 1 - minimal output
 * 2 - intermediate information
 * 3 - very verbose
 * 4 - EVERYTHING
 * */
#define DEBUG 0


/* Struct for Apollonian cell encoder
 * See https://hal.inria.fr/hal-01571226 for theory.
 * Apollonian cell encoders must know:
 * the size of the message space (msize)
 * the size of the key space (ksize)
 * the LCM of the denominators of the irreducible message probabilities (mu)
 * the numerators of the normalised message probabilities (mnums[])
 *     the size of mnums[] should be msize, violating this is unsafe
 * the largest numerator (x)
 * */
struct encoder {
  NUM_TYPE msize;
  NUM_TYPE ksize;
  NUM_TYPE mu;
  NUM_TYPE *mnums;
  NUM_TYPE x;
};


// Support functions

// Greatest common divisor
NUM_TYPE gcd(NUM_TYPE a, NUM_TYPE b) {
  if (DEBUG > 2) {
    printf("GCD("NUM_FMT","NUM_FMT")=",a,b);
  }
  while (a != b) {
    if (a > b) {
      a = a - b;
    } else {
      b = b - a;
    }
  }
  if (DEBUG > 2) {
    printf(""NUM_FMT"\n",a);
  }
  return a;
}

// Lowest common multiple
NUM_TYPE lcm(NUM_TYPE a, NUM_TYPE b) {
  if (DEBUG > 2) {
    NUM_TYPE temp = gcd(a, b);
    temp = (a * b) / temp;
    printf("LCM("NUM_FMT","NUM_FMT")="NUM_FMT"\n",a,b,temp);
    return temp;
  }
  return (a * b) / gcd(a,b); // Division is safe here due to LCM properties
}

// Any initialisation required
void init() {
  // used for rand()
  // NOTE: not cryptographically secure, fix when local_random_sized below is fixed
  srand(time(NULL)); 
}

// Generate random numbers up to specified number of characters
// Argument cs is number of characters (in local architecture) to generate
NUM_TYPE local_random_sized(NUM_TYPE cs) {
  // For proof of concept, use better random here
  NUM_TYPE i = 0;
  NUM_TYPE res = 0;
  NUM_TYPE m = 1 << CHAR_BIT;
  while (i < cs) {
    res = res << CHAR_BIT;
    // TODO: replace with cryptographically secure random
    res = res + ((NUM_TYPE)(rand()) % m);
    i++;
  }
  if (DEBUG > 2) {
    printf("Random number="NUM_FMT"\n",res);
  }
  return res;
}


// Generate random numbers up to size of NUM_TYPE
// Local to not be confused with any other random function
NUM_TYPE local_random() { 
  return local_random_sized(sizeof(NUM_TYPE));
}

// Destroy an encoder by freeing all the memory
void destroy(struct encoder* enc) {
  free(enc->mnums);
  free(enc);
}

// Builds an encoder function from the rational representations
// Note that the pointer fracs is assumed to be pairs of NUM_TYPE's
// where the first is the numberator n and the second the denominator
// d such that each entry is of the form (n/d).
// Also there are "size" possible messages, so fracs contains 2 * size values.
// NOTE: assumes that fracs are irreducible and non-zero!
struct encoder * build_encoder(NUM_TYPE * fracs, NUM_TYPE size) {
  if (size == 0) {
    // Non-sensical, abort!
    return NULL;
  }
  //struct encoder *res = (struct encoder*)malloc(sizeof(struct encoder));
  struct encoder *res = malloc(sizeof(struct encoder));
  if (res == NULL) {
    // Failed allocation, abort!
    return res;
  }
  NUM_TYPE tmp = size * sizeof(NUM_TYPE);
  res->mnums = malloc(tmp);
  if(res->mnums == NULL) {
    // Failed allocation, clean up and then abort!
    free(res);
    return NULL;
  }
  res->mu = 1;
  res->msize = size;
  res->x = 1;
  NUM_TYPE i = 0;
  // First pass, calculate mu
  while (i < size) {
  /*
    // The following calculation of GCD and update ensure fractions are irreducible
    printf("About to do gcd("NUM_FMT","NUM_FMT")\n", fracs[2*i], fracs[2*i+1]);
    NUM_TYPE g = gcd(fracs[2*i], fracs[2*i+1]);
    printf("Found gcd = "NUM_FMT"\n", g);
    if (g != 1) {
      if (TRUE) { //(DEBUG > 2) {
        printf("Found fraction "NUM_FMT"/"NUM_FMT", reducing by "NUM_FMT"...\n",fracs[2*i], fracs[2*i+1],g);
      }
      fracs[2*i] = fracs[2*i]/g;
      fracs[2*i+1] = fracs[2*i+1]/g;
    }
  */
    // Update ongoing LCM to yield mu
    res->mu = lcm(res->mu,fracs[2*i+1]);
    i++;
  }
  // Second pass, populate values
  i = 0;
  NUM_TYPE ctr = 0;
  while(i < size) {
    NUM_TYPE tmp = res->mu/fracs[2*i+1];  // Division safe by LCM properties
    tmp = fracs[2*i]*tmp;
    if (res->x <= tmp) {
      res->x = tmp;
    }
    ctr += tmp;
    res->mnums[i] = ctr;
    i++;
  }
  // Finish by calculating the maximal key size
  res->ksize = res->mu/res->x;  // assume this floors for uint division!!!!
  while (res->ksize * res->x > res->mu) {
    res->ksize = res->ksize - 1;
  }
  // TODO: make sure the above floors correctly
  
  // Sanity check
  if (ctr != res->mu) {
    printf("Invalid distribution, aborting!\n");
    destroy(res);
    return NULL;
  }
  return res;
}

// Builds an encoder (see above) with key size "k".
// Returns NULL if the key size is not possible
struct encoder * build_encoder_k(NUM_TYPE * fracs, NUM_TYPE size, NUM_TYPE ksize) {
  struct encoder* temp = build_encoder(fracs, size);
  // Succeeded in building an encoder, but the key size doesn't work,
  // free memory and return NULL
  if (temp == NULL) {
    return NULL;
  }
  if (temp->ksize < ksize) {
    destroy(temp);
    return NULL;
  }
  temp->ksize = ksize;
  return temp;
}


// Encode with given encoder, message, and key
NUM_TYPE encode(struct encoder* enc, NUM_TYPE m, NUM_TYPE k) {
  // error cases
  // TODO: find a better way to report errors, pass in reference to returned value and return a code?
  if (enc == NULL || m >= enc->msize || k >= enc->ksize) {
    return 0;
  }
  NUM_TYPE b = local_random();
  NUM_TYPE bplus = 0;
  NUM_TYPE range;
  // SIDECHANNEL: below may be vulnerable to cache sidechannel attacks
  if (m == 0) {
    range = enc->mnums[0];
  } else {
    range = enc->mnums[m] - enc->mnums[m-1];
    bplus = enc->mnums[m-1];
  }
  // END_SIDECHANNEL
  b = b % range;  // maybe replace modulo here?
  if (DEBUG > 2) {
    printf("Randomised b is "NUM_FMT"\n",b);
    NUM_TYPE tt = bplus + range -1;
    printf("Size of b is "NUM_FMT" starting from "NUM_FMT" (i.e. "NUM_FMT" to "NUM_FMT")\n",range, bplus, bplus, tt);
  }
  NUM_TYPE res = bplus + b + k * enc->x;
  if (DEBUG > 2) {
    printf("(bplus "NUM_FMT") + (b "NUM_FMT") + (k "NUM_FMT") * (enc->x "NUM_FMT") = "NUM_FMT"\n",bplus,b,k,enc->x, res);
  }
  if (res >= enc->mu) {
    res -= enc->mu;
  }
  return res;
}

// Decode with given encoder, ciphertext, and key
NUM_TYPE decode(struct encoder* enc, NUM_TYPE c, NUM_TYPE k) {
  // error cases
  if (enc == NULL || c >= enc->mu || k >= enc->ksize) {
    return 0;
  }
  if (DEBUG > 3) {
    printf("Decoding c="NUM_FMT" and k="NUM_FMT"...",c,k);
  }
  // Below code does the reverse of the modulo operation
  // Done this way to use unsigned types
  NUM_TYPE km = k * enc->x;
  if (km > c) {
    c += enc->mu;
  }
  NUM_TYPE b = c - km;
  if (DEBUG > 3) {
    printf(" b="NUM_FMT"...",b);
  }
  // Find the right message
  NUM_TYPE i = 0;
  // SIDECHANNEL: below may be vulnerable to timing attacks
  while(TRUE) {
    if (DEBUG >4) {
      printf("\n  testing ("NUM_FMT" < "NUM_FMT")...",b,enc->mnums[i]);
    }
    if (b < enc->mnums[i]) {
      if (DEBUG >3) { 
        printf(" m="NUM_FMT"\n",i);
      }
      return i;
    }
    i++;
  }
  // END_SIDECHANNEL
}

// Printing an encoder
void print_encoder(struct encoder *enc) {
  printf("msize =\t"NUM_FMT"\n",enc->msize);
  printf("ksize =\t"NUM_FMT"\n",enc->ksize);
  printf("mu =\t"NUM_FMT"\n",enc->mu);
  printf("x =\t"NUM_FMT"\n",enc->x);
  NUM_TYPE i = 0;
  NUM_TYPE ctr = 0;
  while (i < enc->msize) {
    printf("\tm"NUM_FMT" b range =\t"NUM_FMT"-"NUM_FMT"\n",i,ctr,(enc->mnums[i])-1);
    ctr = enc->mnums[i];
    i++;
  }
}



// Testing code below

// Number of iterations of the test to do
#define ITERATIONS 100


// Do a full test of the message x ker pairs
void full_test(struct encoder *enc) {
  NUM_TYPE i = 0;
  while (i < enc->msize) {
    NUM_TYPE j = 0;
    while (j < enc->ksize) {
      NUM_TYPE c = encode(enc,i,j);
      NUM_TYPE r = decode(enc,c,j);
      if (r != i) {
        printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and c = "NUM_FMT" and dec("NUM_FMT","NUM_FMT") = "NUM_FMT"\n",i,j,c,c,j,r);
      }
      j++;
    }
    i++;
  }
}

// Test all possible message and key pairs ITERATIONS times
void exhaustive_test(struct encoder *enc) {
  NUM_TYPE i = 0;
  while (i < enc->msize) {
    NUM_TYPE j = 0;
    while (j < enc->ksize) {
      NUM_TYPE x = 0;
      while (x < ITERATIONS) {
        NUM_TYPE c = encode(enc,i,j);
        NUM_TYPE r = decode(enc,c,j);
        if (r != i) {
          printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and c = "NUM_FMT" and dec("NUM_FMT","NUM_FMT") = "NUM_FMT"\n",i,j,c,c,j,r);
        }
        x++;
      }
      j++;
    }
    i++;
  }
}

// A simple distribution
NUM_TYPE dist1[] = {3,10,3,10,1,5,2,10};


void test1() {
  printf("Executing test 1...\n");
  struct encoder* enc = build_encoder_k(dist1,4,3);
  print_encoder(enc);
  NUM_TYPE m;
  NUM_TYPE k;
  NUM_TYPE c;
  NUM_TYPE i;
  // m = 0 tests
  i = 0;
  m = 0;
  k = 0;
  while (i < ITERATIONS) {
    c = encode(enc,m,k);
    if (c > 2) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and result was c = "NUM_FMT"\n",m,k,c);
    }
    i++;
  }
  i = 0;
  m = 0;
  k = 1;
  while (i < ITERATIONS) {
    c = encode(enc,m,k);
    if (c < 3 || c > 5) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and result was c = "NUM_FMT"\n",m,k,c);
    }
    i++;
  }
  i = 0;
  m = 0;
  k = 2;
  while (i < ITERATIONS) {
    c = encode(enc,m,k);
    if (c < 6 || c > 8) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and result was c = "NUM_FMT"\n",m,k,c);
    }
    i++;
  }
  // m = 1 tests
  i = 0;
  m = 1;
  k = 0;
  while (i < ITERATIONS) {
    c = encode(enc,m,k);
    if (c < 3 || c > 5) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and result was c = "NUM_FMT"\n",m,k,c);
    }
    i++;
  }
  i = 0;
  m = 1;
  k = 1;
  while (i < ITERATIONS) {
    c = encode(enc,m,k);
    if (c < 6 || c > 8) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and result was c = "NUM_FMT"\n",m,k,c);
    }
    i++;
  }
  i = 0;
  m = 1;
  k = 2;
  while (i < ITERATIONS) {
    c = encode(enc,m,k);
    if (!(c == 9 || c == 0 || c == 1)) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and result was c = "NUM_FMT"\n",m,k,c);
    }
    i++;
  }
  // m = 2 tests
  i = 0;
  m = 2;
  k = 0;
  while (i < ITERATIONS) {
    c = encode(enc,m,k);
    if (c < 6 || c > 7) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and result was c = "NUM_FMT"\n",m,k,c);
    }
    i++;
  }
  i = 0;
  m = 2;
  k = 1;
  while (i < ITERATIONS) {
    c = encode(enc,m,k);
    if (!(c == 9 || c == 0 )) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and result was c = "NUM_FMT"\n",m,k,c);
    }
    i++;
  }
  i = 0;
  m = 2;
  k = 2;
  while (i < ITERATIONS) {
    c = encode(enc,m,k);
    if (c < 2 || c > 3) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and result was c = "NUM_FMT"\n",m,k,c);
    }
    i++;
  }
  // m = 3 tests
  i = 0;
  m = 3;
  k = 0;
  while (i < ITERATIONS) {
    c = encode(enc,m,k);
    if (c < 8 || c > 9) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and result was c = "NUM_FMT"\n",m,k,c);
    }
    i++;
  }
  i = 0;
  m = 3;
  k = 1;
  while (i < ITERATIONS) {
    c = encode(enc,m,k);
    if (c < 1 || c > 2) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and result was c = "NUM_FMT"\n",m,k,c);
    }
    i++;
  }
  i = 0;
  m = 3;
  k = 2;
  while (i < ITERATIONS) {
    c = encode(enc,m,k);
    if (c < 4 || c > 5) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and result was c = "NUM_FMT"\n",m,k,c);
    }
    i++;
  }
  
  // Random tests
  i = 0;
  while (i < ITERATIONS) {
    m = local_random() % 4;
    k = local_random() % 3;
    c = encode(enc,m,k);
    NUM_TYPE r = decode(enc,c,k);
    if (m != r) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and c = "NUM_FMT" and dec("NUM_FMT","NUM_FMT") = "NUM_FMT"\n",m,k,c,c,k,r);
    }
    i++;
  }

  // clean up
  destroy(enc);
  printf("Finished executing test 1.\n");
}


// The code below does random testing

void print_dist(NUM_TYPE size, NUM_TYPE* vals) {
  NUM_TYPE i = 0;
  while (i < size) {
    printf(""NUM_FMT"/"NUM_FMT"\n", vals[2*i],vals[2*i+1]);
    i++;
  }
}

// Build a not so random distribution
NUM_TYPE* build_dist(NUM_TYPE size) {
  NUM_TYPE* vals;
  if (size < 2) {
    // Non-sensical, do nothing
    return NULL;
  }
  vals = malloc(2 * size * sizeof(NUM_TYPE));
  if (vals == NULL) {
    return NULL;
  }
  // Magic numbers for testing!!!
  NUM_TYPE mu_max = 100;  // soft max, will be pushed up by large size
  // End magic numbers
  NUM_TYPE mu = local_random() % mu_max;
  if (mu <= size) {
    mu = mu + size;
  }
  NUM_TYPE rem = mu;
  NUM_TYPE i = size;
  while (i > 1) {
    i = i - 1;
    if ((rem - i)/2 == 0) {
      vals[2*i] = 1;
    } else {
      vals[2*i] = (local_random() % ((rem - i)/2)) + 1;
    }
    rem = rem - vals[2*i];
    vals[2*i+1] = mu;
  }
  vals[0] = rem;
  vals[1] = mu;
  if (DEBUG > 2) {
    printf("Generated the following distribution:\n");
    print_dist(size,vals);
  }
  return vals;
}

// Run a single random test
void random_test() {
  // Magic testing numbers
  NUM_TYPE max_dsize = 10;
  NUM_TYPE min_dsize = 5;
  // End magic numbers
  NUM_TYPE size = local_random() % max_dsize;
  if (size < min_dsize) {
    size = min_dsize;
  }
  NUM_TYPE *vals;
  // buildDist(size,vals);  // unsafe
  vals = build_dist(size);
  if (vals == NULL) {
    printf("Allocation failed, aborting random test!\n");
  }
  struct encoder* enc = build_encoder(vals,size);
  free(vals);
  if (!(enc)) {
    return;
  }
  if (enc->ksize < 1) {
    printf("Created encoder with valid key size "NUM_FMT", aborting.\n", enc->ksize);
    destroy(enc);
    return;
  }
  if (DEBUG > 0) {
    print_encoder(enc);
  }
  NUM_TYPE i = 0;
  while (i < ITERATIONS) {
    NUM_TYPE m = local_random() % enc->msize;
    NUM_TYPE k = local_random() % enc->ksize;
    NUM_TYPE c = encode(enc,m,k);
    NUM_TYPE r = decode(enc,c,k);
    if (m != r) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and c = "NUM_FMT" and dec("NUM_FMT","NUM_FMT") = "NUM_FMT"\n",m,k,c,c,k,r);
    }
    i++;
  }
  destroy(enc);
}

// Run a single random test
void exhaustive_random_test() {
  // Magic testing numbers
  NUM_TYPE max_dsize = 10;
  NUM_TYPE min_dsize = 5;
  // End magic numbers
  NUM_TYPE size = local_random() % max_dsize;
  if (size < min_dsize) {
    size = min_dsize;
  }
  NUM_TYPE *vals;
  // buildDist(size,vals);  // unsafe
  vals = build_dist(size);
  if (vals == NULL) {
    printf("Allocation failed, aborting random test!\n");
  }
  struct encoder* enc = build_encoder(vals,size);
  free(vals);
  if (!(enc)) {
    return;
  }
  if (enc->ksize < 1) {
    printf("Created encoder with valid key size "NUM_FMT", aborting.\n", enc->ksize);
    destroy(enc);
    return;
  }
  if (DEBUG > 0) {
    print_encoder(enc);
  }
  exhaustive_test(enc);
  destroy(enc);
}

// Run a single random test and show distribution information
void distribution_random_test(NUM_TYPE iters) {
  // Magic testing numbers
  NUM_TYPE max_dsize = 15;
  NUM_TYPE min_dsize = 5;
  // End magic numbers
  NUM_TYPE size = local_random() % max_dsize;
  if (size < min_dsize) {
    size = min_dsize;
  }
  NUM_TYPE *vals;
  // buildDist(size,vals);  // unsafe
  vals = build_dist(size);
  if (vals == NULL) {
    printf("Allocation failed, aborting random test!\n");
    return;
  }
  struct encoder* enc = build_encoder(vals,size);
  if (!(enc)) {
    free(vals);
    return;
  }
  if (enc->ksize < 1) {
    printf("Created encoder with valid key size "NUM_FMT", aborting.\n", enc->ksize);
    destroy(enc);
    return;
  }
  if (DEBUG > 0) {
    print_encoder(enc);
  }
  NUM_TYPE *mc = malloc(enc->msize * sizeof(NUM_TYPE));
  if (mc == NULL) {
    free(vals);
    destroy(enc);
    return;
  }
  NUM_TYPE *cc = malloc(enc->mu * sizeof(NUM_TYPE));
  if (cc == NULL) {
    free(mc);
    free(vals);
    destroy(enc);
    return;
  }
  NUM_TYPE i = 0;
  while (i < enc->msize) {
    mc[i] = 0;
    i++;
  }
  i = 0;
  while (i < enc->mu) {
    cc[i] = 0;
    i++;
  }
  i = 0;
  while (i < iters) {
    NUM_TYPE m = decode(enc,local_random() % enc->mu,0); // hack
    NUM_TYPE k = local_random() % enc->ksize;
    mc[m]++;
    NUM_TYPE c = encode(enc,m,k);
    cc[c]++;
    NUM_TYPE r = decode(enc,c,k);
    if (r != m) {
      printf("ERROR! Had m="NUM_FMT" and k="NUM_FMT" and c = "NUM_FMT" and dec("NUM_FMT","NUM_FMT") = "NUM_FMT"\n",m,k,c,c,k,r);
    }
    i++;
  }
  i = 0;
  printf("Results of message generation distribution:\n");
  while (i < enc->msize) {
    printf(""NUM_FMT": Defined "NUM_FMT"/"NUM_FMT", actual "NUM_FMT"/"NUM_FMT"\n",i,vals[2*i],vals[2*i+1],mc[i],iters);
    i++;
  }
  i = 0;
  printf("Results of ciphertext distribution:\n");
  while (i < enc->mu) {
    printf(""NUM_FMT": Actual "NUM_FMT"/"NUM_FMT"\n",i,cc[i],iters);
    i++;
  }
  free(mc);
  free(cc);
  free(vals);
  destroy(enc);
}

void multi_test(NUM_TYPE count) {
  printf("Performing multi_test for "NUM_FMT" encoders.\n", count);
  NUM_TYPE i = 0;
  while (i < count) {
    i++;
    if (DEBUG > 0 ) {
      printf("Performing random test number "NUM_FMT"... ",i);
    }
    random_test();
    if (DEBUG > 0 ) {
      printf(" finished!\n");
    }
  }
}

void exhaustive_multi_test(NUM_TYPE count) {
  printf("Performing exhaustive_multi_test for "NUM_FMT" encoders.\n", count);
  NUM_TYPE i = 0;
  while (i < count) {
    i++;
    if (DEBUG > 0 ) {
      printf("Performing exhaustive random test number "NUM_FMT"... ",i);
    }
    exhaustive_random_test();
    if (DEBUG > 0 ) {
      printf(" finished!\n");
    }
  }
}

void distribution_multi_test(NUM_TYPE count, NUM_TYPE iters) {
  printf("Performing distribution_multi_test for "NUM_FMT" encoders.\n", count);
  NUM_TYPE i = 0;
  while (i < count) {
    i++;
    if (DEBUG > 0 ) {
      printf("Performing distribution random test number "NUM_FMT"... ",i);
    }
    distribution_random_test(iters);
    if (DEBUG > 0 ) {
      printf(" finished!\n");
    }
  }
}

/* CISCO DEMO CODE */
void run_demo() {
  NUM_TYPE demo_dist[] = {3,10,3,10,1,5,2,10};
  NUM_TYPE demo_hack[] = {1,4,1,4,1,4,1,4};
  printf("Building Latin Encoder...\n");
  struct encoder* la_enc = build_encoder_k(demo_hack,4,3);
  printf("Building Apollonian Cell Encoder...\n");
  struct encoder* ap_enc = build_encoder_k(demo_dist,4,3);
  NUM_TYPE k = local_random() % 3;
  printf("Choosing key randomly; chose "NUM_FMT"\n",k);
  NUM_TYPE i = 0;
  NUM_TYPE ms[4];
  NUM_TYPE cs[10];
  while (i < 4) {
    ms[i] = 0;
    i++;
  }
  i = 0;
  while (i < 10) {
    cs[i] = 0;
    i++;
  }
  i = 0;
  printf("Simulating encoding 10000 messages with Latin encoder.\n");
  while (i < 10000) {
    NUM_TYPE m = local_random() % 10;
    if (m < 3) {
      m = 0;
    } else if (m < 6) {
      m = 1;
    } else if (m < 8) {
      m = 2;
    } else {
      m = 3;
    }
    ms[m]++;
    NUM_TYPE c = encode(la_enc,m,k);
    cs[c]++;
    // Could validate here ...
    i++;
  }
  i = 0;
  printf("Actual message counts produced:\n");
  while (i < 4) {
    printf("Message "NUM_FMT" = "NUM_FMT"\n",i,ms[i]);
    i++;
  }
  i = 0;
  printf("Actual ciphertext counts produced:\n");
  while (i < 4) {
    printf("Ciphertext "NUM_FMT" = "NUM_FMT"\n",i,cs[i]);
    i++;
  }
  printf("Making guess about key... (note numbers are relative)\n");
  // Insert proper calculations here
  NUM_TYPE ek0 = 24* cs[0] + 21*cs[1] + 14* cs[2] + 16*cs[3];
  NUM_TYPE ek1 = 16* cs[0] + 21*cs[1] + 21* cs[2] + 16*cs[3];
  NUM_TYPE ek2 = 16* cs[0] + 14*cs[1] + 21* cs[2] + 24*cs[3];
  NUM_TYPE kg = 0;
  if (ek0 > ek1) {
    printf("\tevidence k0 ("NUM_FMT") > k1 ("NUM_FMT")...\n", ek0, ek1);
    if (ek0 > ek2) {
      printf("\tevidence k0 ("NUM_FMT") > k2 ("NUM_FMT")...\nResult is k=0\n", ek0, ek2);
      kg = 0;
    } else {
      printf("\tevidence k0 ("NUM_FMT") <= k2 ("NUM_FMT")...\nResult is k=2\n", ek0, ek2);
      kg = 2;
    }
  } else {
    printf("\tevidence k0 ("NUM_FMT") <= k1 ("NUM_FMT")...\n", ek0, ek1);
    if (ek1 > ek2) {
      printf("\tevidence k1 ("NUM_FMT") > k2 ("NUM_FMT")...\nResult is k=1\n", ek1, ek2);
      kg = 1;
    } else {
      printf("\tevidence k1 ("NUM_FMT") <= k2 ("NUM_FMT")...\nResult is k=2\n", ek1, ek2);
      kg = 2;
    }
  }
  if (kg == k) {
    printf("CORRECT!\n");
  } else {
    printf("INCORRECT!\n");
  }
  /* Old approximation below
  if (cs[1] > cs[2] + 500) {
    printf("since "NUM_FMT" >> "NUM_FMT" guessing key is 0.\n", cs[1],cs[2]);
  } else if (cs[2] > cs[3] + 500) {
    printf("since "NUM_FMT" >> "NUM_FMT" guessing key is 1.\n", cs[2],cs[3]);
  } else {
    printf("since "NUM_FMT" >> "NUM_FMT" guessing key is 2.\n", cs[3],cs[0]);
  }*/
  i = 0;
  while (i < 4) {
    ms[i] = 0;
    i++;
  }
  i = 0;
  while (i < 10) {
    cs[i] = 0;
    i++;
  }
  i = 0;
  printf("Simulating encoding 10000 messages with Apollonian Cell Encoder.\n");
  while (i < 10000) {
    NUM_TYPE m = local_random() % 10;
    if (m < 3) {
      m = 0;
    } else if (m < 6) {
      m = 1;
    } else if (m < 8) {
      m = 2;
    } else {
      m = 3;
    }
    ms[m]++;
    NUM_TYPE c = encode(ap_enc,m,k);
    cs[c]++;
    // Could validate here ...
    i++;
  }
  i = 0;
  printf("Actual message counts produced:\n");
  while (i < 4) {
    printf("Message "NUM_FMT" = "NUM_FMT"\n",i,ms[i]);
    i++;
  }
  i = 0;
  printf("Actual ciphertext counts produced:\n");
  while (i < 10) {
    printf("Ciphertext "NUM_FMT" = "NUM_FMT"\n",i,cs[i]);
    i++;
  }
  printf("Making guess about key... (note numbers are relative)\n");
  // insert proper calculations here
  ek0 = 0;
  i = 0;
  while (i < 10) {
    ek0 += cs[i];
    i++;
  }
  ek1 = ek0;
  ek2 = ek0;
  kg = local_random() %3;
  printf("Guessing "NUM_FMT" since evidence for k=0 ("NUM_FMT") == k=1 ("NUM_FMT") == k=2 ("NUM_FMT").\n",kg,ek0,ek1,ek2);
  if (kg == k) {
    printf("CORRECT!\n");
  } else {
    printf("INCORRECT!\n");
  }
  printf("Finished simulation, cleaning up.\n");
  destroy(la_enc);
  destroy(ap_enc);
}

/* END DEMO CODE */

int main() {
  init();
  test1();
  // randomTest(); 
  multi_test(1000);
  exhaustive_multi_test(1000);
  distribution_multi_test(10,100000);
  //run_demo();
  exit(0);
}


