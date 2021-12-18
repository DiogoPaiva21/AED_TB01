//
// AED, November 2021
//
// Solution of the first practical assignement (subset sum problem)
//
// Place your student numbers and names here
//

#if __STDC_VERSION__ < 199901L
#error "This code must must be compiled in c99 mode or later (-std=c99)" // to handle the unsigned long long data type
#endif
#ifndef STUDENT_H_FILE
#define STUDENT_H_FILE "000000.h"
#endif

//
// include files
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "elapsed_time.h"
#include STUDENT_H_FILE

//
// custom data types
//
// the STUDENT_H_FILE defines the following constants and data types
//
//   #define min_n       24                   --- the smallest n value we will handle
//   #define max_n       57                   --- the largest n value we will handle
//   #define n_sums      20                   --- the number of sums for each n value
//   #define n_problems  (max_n - min_n + 1)  --- the number of n values
//
//   typedef unsigned long long integer_t;    ---  64-bit unsigned integer
//   typedef struct
//   {
//     int n;                                 --- number of elements of the set (for a valid problem, min_n <= n <= max_n)
//     integer_t p[max_n];                    --- the elements of the set, already sorted in increasing order (only the first n elements are used)
//     integer_t sums[n_sums];                --- several sums (problem: for each sum find the corresponding subset)
//   }
//   subset_sum_problem_data_t;               --- weights p[] and sums for a given value of n
//
//   subset_sum_problem_data_t all_subset_sum_problems[n_problems]; --- // the problems
//

//
// place your code here
//
// possible function prototype for a recursive brute-force function:
// int brute_force(int n,integer_t p[n],int desired_sum,int current_index,integer_t partial_sum)
// it sould return 1 when the solution is found and 0 otherwise
// note, however, that you may get a faster function by reducing the number of function arguments (maybe a single pointer to a struct?)
//

/**
 * @brief Determina a combinação dos valores de p cujo somatório é desired_sum, por um método brute force não recursivo.
 *
 * @param n             Tamanho de p
 * @param p             Conjunto com os valores a somar
 * @param desired_sum   Soma desejada
 * @param b             Array que guarda a solução
 *
 * @return              1 se foi encontrada uma solução, 0 caso contrário
 */
int brute_force(int n, integer_t p[n], integer_t desired_sum, int b[n])
{
    integer_t mask;     // máscara
    integer_t test_sum; // soma de teste

    /* para cada combinação */
    for (mask = 0; mask < (1 << n); mask++)
	{
        /* determinar a soma dos valores de p */
        test_sum = 0;
        for (int i = 0; i < n; i++)
            if ((mask & (1 << i)) != 0)
                test_sum += p[i];

        /* foi descoberta a soma */
		if (test_sum == desired_sum)
		{
            /* guardar a máscara no array b */
            for (int i = 0; i < n; i++)
                b[i] = ((mask & (1 << i)) == 0) ? 0 : 1;
			return 1;
		}
	}

	return 0;
}

/**
 * @brief Determina a combinação dos valores de p cujo somatório é desired_sum, por um método brute force recursivo.
 *
 * @param n                 Tamanho de p
 * @param p                 Conjunto com os valores a somar
 * @param desired_sum       Soma desejada
 * @param current_index     Índice atual em p
 * @param partial_sum       Soma atual
 * @param mask              Máscara atual
 * @param b                 Array que guarda a solução
 *
 * @return                  1 se foi encontrada uma solução, 0 caso contrário
 *
 */
int brute_force_recursive(int n, integer_t p[n], integer_t desired_sum, int current_index, integer_t partial_sum, integer_t mask, int b[n])
{
    /* foi descoberta a soma */
	if (partial_sum == desired_sum)
	{
        /* guardar a máscara no array b */
        for (int i = 0; i < n; i++)
            b[i] = ((mask & (1 << i)) == 0) ? 0 : 1;
		return 1;
	}

    /* valor de retorno é zero se ainda não tiver sido encontrada uma solução */
    int ret = 0;
    if (current_index < n)
    {
        /* p[current_index] não é usado na soma */
        ret |= brute_force_recursive(n, p, desired_sum, current_index + 1, partial_sum, mask, b);
        /* p[current_index] é usado na soma */
        ret |= brute_force_recursive(n, p, desired_sum, current_index + 1, partial_sum + p[current_index], mask | (1 << current_index), b);
	}
    return ret;
}

//
// main program
//

int main(void)
{
	fprintf(stderr, "Program configuration:\n");
	fprintf(stderr, "  min_n ....... %d\n", min_n);
	fprintf(stderr, "  max_n ....... %d\n", max_n);
	fprintf(stderr, "  n_sums ...... %d\n", n_sums);
	fprintf(stderr, "  n_problems .. %d\n", n_problems);
	fprintf(stderr, "  integer_t ... %d bits\n", 8 * (int)sizeof(integer_t));

	//
    // para cada problema
    //
	for (int i = 0; i < n_problems; i++)
	{
		int n = all_subset_sum_problems[i].n;           // numero de valores a somar		
		integer_t *p = all_subset_sum_problems[i].p;    // valores a somar

		if (n > 10) continue;   // ignorar problemas com numero de valores a somar superior

		//
        // para cada soma
        //
		for (int j = 0; j < n_sums; j++)
		{
			integer_t desired_sum = all_subset_sum_problems[i].sums[j]; // soma desejada
			int b[n];   // resultados

            /* 
             * brute force não recursiva
             */
			//brute_force(all_subset_sum_problems[i].n, all_subset_sum_problems[i].p, all_subset_sum_problems[i].sums[j]);

            /* 
             * brute force recursiva
             */
			brute_force_recursive(n, p, desired_sum, 0, 0, b);

			// imprime o resultado
			for (int i = 0; i < n; i++)
				printf("%d", b[i]);
			printf("\n");
		}
	}

	return 0;
}
