﻿//Безпалый Денис КИ22-06Б
//Вариант 7
//Практическая работа 1

#include <iostream>
#include <cmath>
#include "omp.h"

using namespace std;

void printArray(double arr[], int size) {
	cout << "Массив: \n";
	for (size_t i = 0; i < size; i++)
	{
		cout << arr[i] << " ";
	}
	cout << "\n";
}

void sequentalCompute(double a[], double b[], double c[], double d[], int N, int x, int L)
{
	a[0] = x;
	b[0] = x;
	d[0] = x;

	//Заполнение массивов
	for (int j = 0; j < L; j++)
	{
		for (int i = 1; i < N; i++) {
			a[i] = a[i - 1] * i;
		}
	}
	for (int j = 0; j < L; j++) {
		for (int i = 1; i < N; i++)
		{
			b[i] = cos(x / i) + x / i;
		}
	}
	for (int j = 0; j < L; j++) {
		for (int i = 0; i < N; i++)
		{
			c[i] = a[i] - b[N - i];
		}
	}
	for (int j = 0; j < L; j++) {
		for (int i = 1; i < N; i++)
		{
			d[i] = sin(d[i - 1] * i);
		}
	}

	//Переопределение массивов
	for (int j = 0; j < L; j++) {
		for (int i = 0; i < N; i++)
		{
			b[i] = (a[i] + c[i]) / 2;
		}
	}
	d[0] = x;
	for (int j = 0; j < L; j++) {
		for (int i = 1; i < N; i++)
		{
			d[i] = sin(a[i] * i) + d[i] / i;
		}
	}
}

void parallelCompute(double a[], double b[], double c[], double d[], int N, int x, int k, int L)
{
	a[0] = x;
	b[0] = x;
	d[0] = x;

	omp_set_num_threads(k);

	// Вычисление массива a
		/*for (int i = 1; i < N; i++) {
			for (int j = 0; j < L; j++) {
				a[i] = a[i - 1] * i;
			}
		}*/

#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
		int chunk = (N - 1) / k + ((N - 1) % k != 0);
		int start_index = thread_id * chunk + 1;
		int end_index = min(start_index + chunk, N);

		for (int i = start_index; i < end_index; i++) {
			for (int j = 0; j < L; j++) {
				a[i] = a[0] * tgamma(i + 1);
			}
		}
	}

	// Вычисление массива b (первый раз)
#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
		int chunk = (N - 1) / k + ((N - 1) % k != 0);
		int start_index = thread_id * chunk + 1;
		int end_index = min(start_index + chunk, N);

		for (int i = start_index; i < end_index; i++) {
			for (int j = 0; j < L; j++) {
				b[i] = cos(x / i) + x / i;
			}
		}
	}

	// Вычисление массива c
#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
		int chunk = N / k + (N % k != 0);
		int start_index = thread_id * chunk;
		int end_index = min(start_index + chunk, N);

		for (int i = 0; i < end_index; i++) {
			for (int j = 0; j < L; j++) {
				c[i] = a[i] - b[N - i];
			}
		}
	}

	// Вычисление массива d (первый раз)
		for (int i = 1; i < N; i++) {
			for (int j = 0; j < L; j++) {
				d[i] = sin(d[i - 1] * i);
			}
		}

	// Переопределение массива b
#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
		int chunk = N / k + (N % k != 0);
		int start_index = thread_id * chunk;
		int end_index = min(start_index + chunk, N);

		for (int i = start_index; i < end_index; i++) {
			for (int j = 0; j < L; j++) {
				b[i] = (a[i] + c[i]) / 2;
			}
		}
	}

	// Переопределение массива d
	d[0] = x;
#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
		int chunk = (N - 1) / k + ((N - 1) % k != 0);
		int start_index = thread_id * chunk + 1;
		int end_index = min(start_index + chunk, N);

		for (int i = start_index; i < end_index; i++) {
			for (int j = 0; j < L; j++) {
				d[i] = sin(a[i] * i) + d[i] / i;
			}
		}
	}
}

int main()
{
	setlocale(0, "");
	int N = 10, k, L = 100;
	double x = 5;
	/*cout << "Введите N: ";
	cin >> N;*/
	cout << "Введите k: ";
	cin >> k;
	//cout << "Введите x: ";
	//cin >> x;

	double* a = new double[N];
	double* b = new double[N];
	double* c = new double[N];
	double* d = new double[N];

	double start, end, time;

	start = omp_get_wtime();
	sequentalCompute(a, b, c, d, N, x, L);
	end = omp_get_wtime();

	time = end - start;

	cout << "Время последовательного вычисления: " << fixed << time << " сек." << endl;

	printArray(a, N);
	printArray(b, N);
	printArray(c, N);
	printArray(d, N);

	start = omp_get_wtime();
	parallelCompute(a, b, c, d, N, x, k, L);
	end = omp_get_wtime();

	time = end - start;

	cout << "Время параллельного вычисления: " << fixed << time << " сек." << endl;

	printArray(a, N);
	printArray(b, N);
	printArray(c, N);
	printArray(d, N);

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] d;

	return 0;
}