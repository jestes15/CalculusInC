#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif

#define __DELTA__ 0.00001

struct threadpool_data {
    double upper_bound;
    double lower_bound;

    double result;

    double (*integrand)(double);
};

struct threadpool_private_data {
    double _thread_one_results;
    double _thread_two_results;
    double _thread_three_results;
    double _thread_four_results;

    double _first_term;
    double _f_x0;
    double _f_xN;

};

struct thread_private_data
{
    double upper_bound;
    double lower_bound;

    double result;

    double (*integrand)(double);
};

void* integral(void *);
double func(double);
void schedule(struct threadpool_data*);
void *get_aux_data(struct threadpool_private_data *, double (*)(double), double, double);

    int main(void)
{
    struct threadpool_data *t1 = (struct threadpool_data *)malloc(sizeof(struct threadpool_data));
    t1->upper_bound = 16;
    t1->lower_bound = 0;
    t1->integrand = func;
    t1->result = 0;

    schedule(t1);

    printf("%g", t1->result);

    free(t1);
}

double func(double x) {
    return x;
}

void schedule(struct threadpool_data *data)
{
    struct thread_private_data* private_data[NUM_THREADS];

    for (int i = 0; i < NUM_THREADS; ++i) {
        private_data[i] = (struct thread_private_data *)malloc(sizeof(struct thread_private_data));
    }

    // Create struct to pass to get_aux_data
    struct threadpool_private_data *aux_data = (struct threadpool_private_data *)malloc(sizeof(struct threadpool_private_data));

    // Get offset for seperating data
    double offset = (data->upper_bound - data->lower_bound) / NUM_THREADS;

    private_data[0]->lower_bound = data->lower_bound + __DELTA__;
    private_data[0]->upper_bound = private_data[0]->lower_bound + offset;

    for (int i = 1; i < NUM_THREADS - 1; ++i) {
        private_data[i]->lower_bound = private_data[i - 1]->upper_bound + __DELTA__;
        private_data[i]->upper_bound = private_data[i]->lower_bound + offset;
    }

    private_data[NUM_THREADS - 1]->lower_bound = private_data[NUM_THREADS - 2]->upper_bound + __DELTA__;
    private_data[NUM_THREADS - 1]->upper_bound = private_data[NUM_THREADS - 1]->lower_bound + offset - __DELTA__;

    pthread_t threads[NUM_THREADS];
    for (int i = 0; i < NUM_THREADS; ++i)
        pthread_create(&threads[i], NULL, integral, (void*)private_data[i]);

    get_aux_data(aux_data, data->integrand, data->upper_bound, data->lower_bound);

    int return_values[NUM_THREADS];
    double _f_sum = 0;

    for (int i = 0; i < NUM_THREADS; ++i) {
        pthread_join(threads[i], (void*)&return_values[i]);
        _f_sum += private_data[i]->result;
    }

    data->result = aux_data->_first_term * (aux_data->_f_x0 + _f_sum + aux_data->_f_xN);

    // Free memory
    free(aux_data);
}

void* get_aux_data(struct threadpool_private_data* data, double (*integrand)(double), double upper_bound, double lower_bound)
{
    data->_first_term = __DELTA__ / 2;
    data->_f_x0 = integrand(lower_bound);
    data->_f_xN = integrand(upper_bound);

    return (void*)0;
}

/**
 * @param __args Is a pointer to the struct threadpool_data
 */
void* integral(void* _T_args)
{
    struct thread_private_data *data = (struct thread_private_data *)_T_args;
    double (*integrand)(double) = data->integrand;

    double _f_sum;

    _f_sum = 0;
    for (double i = data->lower_bound; i <= data->upper_bound; i += __DELTA__)
        _f_sum += 2 * integrand(i);

    data->result = _f_sum;

    return (void*)EXIT_SUCCESS;
}
