// Location.cpp : ���� DLL Ӧ�ó���ĵ���������
//
#include "stdio.h"
#include "math.h"

#define N 6
#define ERR_SUCCESS	0    //�ɹ�
#define ERR_INPUTA	1    //�����������
#define ERR_INVERSION	2    //����ʧ��
#define ERR_FAIL	3    //ʧ��


int Sanbian(float(*MA)[N], float(*MB)[N], float(*MR)[N]);
void Multiplication(float(*MA)[N], float(*MB)[N], float(*MR)[N]);
int Inversion(float(*MA)[N], float(*MR)[N]);
void Transposition(float(*MA)[N], float(*MR)[N]);
float getA(float(*arcs)[N], int n);
void getAStart(float(*arcs)[N], int n, float(*ans)[N]);

int GetPositionAuto(float(*inputA)[4], float(*outputA), int n);
int GetPosition1D(float(*inputA)[4], float(*outputA), int n);
float GetPosition2D(float(*inputA)[4], float(*outputA), int n);
int GetPosition2DNew(float(*inputA)[4], float(*outputA), int n);
int GetPosition3D(float(*inputA)[4], float(*outputA), int n);
float GetResidualError(float(*inputA)[4], float(*outputA), int n);


int point_in_triangle(float(*input1)[4], float(*input2)[4], int m, int n, int k);
int point_in_rectangle(float(*input1)[4], float(*input2)[4]);
int return_value(float(*input1)[4], float(*input2)[4], float(*output)[4]);
//���淽���ά��
int dimensional = 2;

int GetPosition(float(*inputA)[4], float(*outputA), int n, int c) {
	int err = ERR_FAIL;
	if (n<2 || n > N) return ERR_INPUTA;
	switch (c)
	{
	case 0: err = GetPositionAuto(inputA, outputA, n); break;
	case 1: err = GetPosition1D(inputA, outputA, n); break;
	case 2: err = GetPosition2D(inputA, outputA, n); break;
	case 3: err = GetPosition3D(inputA, outputA, n); break;
	case 4: err = GetPosition2DNew(inputA, outputA, n); break;
	default:
		break;
	}
	return err;
}

//Anchor����Ϊ4����ʱ����3D��3ʱ����2D,2ʱ����1D
int GetPositionAuto(float(*inputA)[4], float(*outputA), int n) {
	int i, j, err;
	float A[N][N] = { 0 }, B[N][N] = { 0 }, C[N][N] = { 0 };
	switch (n)
	{
	case 2:err = GetPosition1D(inputA, outputA, n);
		break;

	case 3:for (i = 0; i < n - 1; i++) {
		for (j = 0; j < 2; j++) {
			//����A����
			*(*(A + i) + j) = 2 * (*(*(inputA + i) + j) - *(*(inputA + n - 1) + j));
		}
		//����B����
		B[i][0] = pow(**(inputA + i), 2) - pow(**(inputA + n - 1), 2) + pow(*(*(inputA + i) + 1), 2) - pow(*(*(inputA + n - 1) + 1), 2)
			+ pow(*(*(inputA + n - 1) + 3), 2) - pow(*(*(inputA + i) + 3), 2);
	}
		   dimensional = 2;
		   //����C����
		   err = Sanbian(A, B, C);

		   for (i = 0; i < 3; i++) {
			   outputA[i] = C[i][0];
		   }
		   break;

	default:for (i = 0; i < n - 1; i++) {
		for (j = 0; j < 3; j++) {
			*(*(A + i) + j) = 2 * (*(*(inputA + i) + j) - *(*(inputA + n - 1) + j));
		}
		B[i][0] = pow(**(inputA + i), 2) - pow(**(inputA + n - 1), 2) + pow(*(*(inputA + i) + 1), 2) - pow(*(*(inputA + n - 1) + 1), 2)
			+ pow(*(*(inputA + i) + 2), 2) - pow(*(*(inputA + n - 1) + 2), 2) + pow(*(*(inputA + n - 1) + 3), 2) - pow(*(*(inputA + i) + 3), 2);
	}
			dimensional = 3;
			//����C����
			err = Sanbian(A, B, C);

			for (i = 0; i < 3; i++) {
				outputA[i] = C[i][0];
			}
			break;
	}
	return err;
}

//һά�㷨
int GetPosition1D(float(*inputA)[4], float(*outputA), int n) {
	float dr = 0, dt = 0, t;
	//�����������ʵ�ʾ���Ͳ��Ծ���
	dr = sqrt(pow(*(*(inputA + 1) + 0) - *(*(inputA + 0) + 0), 2) + pow(*(*(inputA + 1) + 1) - *(*(inputA + 0) + 1), 2) + pow(*(*(inputA + 1) + 2) - *(*(inputA + 0) + 2), 2));
	//Ĭ�ϴ����������֮��
	dt = *(*(inputA + 0) + 3) - (*(*(inputA + 0) + 3) + *(*(inputA + 1) + 3) - dr) * 0.5;
	//�������(��)��
	if (*(*(inputA + 1) + 3) > dr&&*(*(inputA + 0) + 3) < dr)
		dt = 0 - (*(*(inputA + 0) + 3) + *(*(inputA + 1) + 3) - dr) * 0.5;
	//�������(��)��
	if (*(*(inputA + 0) + 3) > dr&&*(*(inputA + 1) + 3) < dr)
		dt = dr + (*(*(inputA + 0) + 3) + *(*(inputA + 1) + 3) - dr) * 0.5;
	//����t
	if (dr != 0) {
		t = dt / dr;
	}
	else {
		return ERR_INPUTA;
	}
	*(outputA + 0) = *(*(inputA + 0) + 0) + t*(*(*(inputA + 1) + 0) - *(*(inputA + 0) + 0));
	*(outputA + 1) = *(*(inputA + 0) + 1) + t*(*(*(inputA + 1) + 1) - *(*(inputA + 0) + 1));
	*(outputA + 2) = *(*(inputA + 0) + 2) + t*(*(*(inputA + 1) + 2) - *(*(inputA + 0) + 2));
	*(outputA + 3) = GetResidualError(inputA, outputA, n);
	return ERR_SUCCESS;
}

//��ά�㷨
float GetPosition2D(float(*inputA)[4], float(*outputA), int n)
{
	int i, j;
	float A[N][N] = { 0 }, B[N][N] = { 0 }, C[N][N] = { 0 };
	if (n < 3) return ERR_INPUTA;
	for (i = 0; i < n - 1; i++) {
		for (j = 0; j < 2; j++) {
			//����A����
			*(*(A + i) + j) = 2 * (*(*(inputA + i) + j) - *(*(inputA + n - 1) + j));
		}
		//����B����
		B[i][0] = pow(**(inputA + i), 2) - pow(**(inputA + n - 1), 2) + pow(*(*(inputA + i) + 1), 2) - pow(*(*(inputA + n - 1) + 1), 2)
			+ pow(*(*(inputA + n - 1) + 3), 2) - pow(*(*(inputA + i) + 3), 2);
	}
	dimensional = 2;
	//����C����
	int err = Sanbian(A, B, C);
	for (i = 0; i < 2; i++) {
		outputA[i] = C[i][0];
	}
	*(outputA + 3) = GetResidualError(inputA, outputA, n);

	return *(outputA + 3);
}

//��ά�㷨
int GetPosition3D(float(*inputA)[4], float(*outputA), int n) {
	int i, j;
	float A[N][N] = { 0 }, B[N][N] = { 0 }, C[N][N] = { 0 };
	if (n < 4) return ERR_INPUTA;
	for (i = 0; i < n - 1; i++) {
		for (j = 0; j < 3; j++) {
			*(*(A + i) + j) = 2 * (*(*(inputA + i) + j) - *(*(inputA + n - 1) + j));
		}
		B[i][0] = pow(**(inputA + i), 2) - pow(**(inputA + n - 1), 2) + pow(*(*(inputA + i) + 1), 2) - pow(*(*(inputA + n - 1) + 1), 2)
			+ pow(*(*(inputA + i) + 2), 2) - pow(*(*(inputA + n - 1) + 2), 2) + pow(*(*(inputA + n - 1) + 3), 2) - pow(*(*(inputA + i) + 3), 2);
	}
	dimensional = 3;
	//����C����
	int err = Sanbian(A, B, C);
	for (i = 0; i < 3; i++) {
		outputA[i] = C[i][0];
	}
	*(outputA + 3) = GetResidualError(inputA, outputA, n);
	return err;
}

int Sanbian(float(*MA)[N], float(*MB)[N], float(*MC)[N]) {
	//C=((A'*A)^(-1))*A'*B
	//MR ���ؽ������ֵ ������Sanbian�㷨��һ����������
	float MR1[N][N] = { 0 };//A'
	float MR2[N][N] = { 0 };//A'*A
	float MR3[N][N] = { 0 };//(A'*A)^(-1)
	float MR4[N][N] = { 0 };//(A'*A)^(-1)*A'

	Transposition(MA, MR1);//A'
	Multiplication(MR1, MA, MR2);//A'*A
	int err = Inversion(MR2, MR3);//(A'*A)^(-1)
	Multiplication(MR3, MR1, MR4);//(A'*A)^(-1)*A'
	Multiplication(MR4, MB, MC);//(A'*A)^(-1)*A'*B

	return err;
}

//�в��㷨
float GetResidualError(float(*inputA)[4], float(*outputA), int n) {
	int i;
	float sum = 0;
	for (i = 0; i < n; i++) {
		sum += pow(sqrt(pow(*(outputA)-*(*(inputA + i)), 2) + pow(*(outputA + 1) - *(*(inputA + i) + 1), 2) + pow(*(outputA + 2) - *(*(inputA + i) + 2), 2)) - *(*(inputA + i) + 3), 2);
	}
	return sqrt(sum / (float)n);
}

//��ά͹�ı����㷨
int GetPosition2DNew(float(*inputA)[4], float(*outputA), int n) {
	int i, j;
	float residuala, residualb, sum1 = 0, sum2 = 0, h = *(outputA + 2);
	float A[N][N] = { 0 }, B[N][N] = { 0 }, C[N][N] = { 0 };
	float Tri[6][4] = { 0 };
	float inputB[1][4] = { 0 }, inputA1[3][4] = { 0 }, inputA2[3][4] = { 0 }, outputA1[4] = { 0 }, outputA2[4] = { 0 };
	if (n < 3) return ERR_INPUTA;
	outputA1[2] = h;
	outputA2[2] = h;
	//��Z��Ͷ����ͬһˮƽ��
	for (i = 0; i < n; i++) {
		*(*(inputA + i) + 3) = sqrt(pow(*(*(inputA + i) + 3), 2) + pow(fabs(*(*(inputA + i) + 2) - h), 2));
	}
	for (i = 0; i < n - 1; i++) {
		for (j = 0; j < 2; j++) {
			//����A����
			*(*(A + i) + j) = 2 * (*(*(inputA + i) + j) - *(*(inputA + n - 1) + j));
		}
		//����B����
		B[i][0] = pow(**(inputA + i), 2) - pow(**(inputA + n - 1), 2) + pow(*(*(inputA + i) + 1), 2) - pow(*(*(inputA + n - 1) + 1), 2)
			+ pow(*(*(inputA + n - 1) + 3), 2) - pow(*(*(inputA + i) + 3), 2);
	}
	dimensional = 2;
	//����C����
	int err = Sanbian(A, B, C);
	for (i = 0; i < 2; i++) {
		outputA[i] = C[i][0];
		inputB[0][i] = C[i][0];
	}
	inputB[0][2] = C[2][0];
	*(outputA + 3) = GetResidualError(inputA, outputA, n);

	if (n == 4) {
		//������������
		return_value(inputA, inputB, Tri);

		for (i = 0; i < 3; i++) {
			for (j = 0; j < 4; j++) {
				sum1 += Tri[i][j];
				*(*(inputA1 + i) + j) = Tri[i][j];
			}
		}
		if (sum1 != 0)
			residuala = GetPosition2D(inputA1, outputA1, 3);

		for (i = 0; i < 3; i++) {
			for (j = 0; j < 4; j++) {
				sum2 += Tri[i + 3][j];
				*(*(inputA2 + i) + j) = Tri[i + 3][j];
			}
		}
		if (sum2 != 0)
			residualb = GetPosition2D(inputA2, outputA2, 3);

		if (sum1 != 0 && sum2 != 0) {
			if (residuala > residualb) {
				for (i = 0; i < 4; i++) {
					outputA[i] = outputA2[i];
				}
			}
			else
			{
				for (i = 0; i < 4; i++) {
					outputA[i] = outputA1[i];
				}
			}
		}
	}
	return err;
}

//������˺���
void Multiplication(float(*MA)[N], float(*MB)[N], float(*MR)[N]) {
	int i, j, k;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				*(*(MR + i) + j) = *(*(MR + i) + j) + *(*(MA + i) + k) * *(*(MB + k) + j);
			}
		}
	}
}

//����ת�ú���
void Transposition(float(*MA)[N], float(*MR)[N]) {
	int i, j;
	for (i = 0; i < N; i++)//������MA�ĵ�n�е�ֵ��������MD�ĵ�n��	
	{
		for (j = 0; j < N; j++)
		{
			MR[j][i] = MA[i][j];
		}
	}
}

//��������溯��
int Inversion(float(*MA)[N], float(*MR)[N]) {
	float arcs[N][N];
	float astar[N][N];
	int i, j;
	//3*3����
	int n = dimensional;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			*(*(arcs + i) + j) = *(*(MA + i) + j);
		}
	}

	float a = getA(arcs, n);
	if (a == 0) {
		return ERR_INVERSION;
		//������ ��������
	}
	else {
		getAStart(arcs, n, astar);

		for (j = 0; j < n; j++) {
			for (i = 0; i < n; i++) {
				*(*(MR + j) + i) = *(*(astar + j) + i) / a;
			}
		}
		return ERR_SUCCESS;
	}
}

//����һ��չ������|A|
float getA(float(*arcs)[N], int n) {
	if (n == 1) {
		return **arcs;
	}
	float ans = 0;
	float temp[N][N];
	int i, j, k;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n - 1; j++) {
			for (k = 0; k < n - 1; k++) {
				*(*(temp + j) + k) = *(*(arcs + j + 1) + ((k >= i) ? k + 1 : k));
			}
		}
		float t = getA(temp, n - 1);
		if (i % 2 == 0) {
			ans += *(*arcs + i) * t;
		}
		else {
			ans -= *(*arcs + i) * t;
		}
	}
	return ans;
}

//����ÿһ��ÿһ�е�ÿ��Ԫ������Ӧ������ʽ�����A*
void getAStart(float(*arcs)[N], int n, float(*ans)[N]) {
	if (n == 1) {
		**ans = 1;
		return;
	}
	int i, j, k, t;
	float temp[N][N];
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n - 1; k++) {
				for (t = 0; t < n - 1; t++) {
					*(*(temp + k) + t) = *(*(arcs + (k >= i ? k + 1 : k)) + (t >= j ? t + 1 : t));
				}
			}
			*(*(ans + j) + i) = getA(temp, n - 1);
			if ((i + j) % 2 == 1) {
				*(*(ans + j) + i) = -*(*(ans + j) + i);
			}
		}
	}
}

//return value:	-1�������������ڣ�0���������α��ϣ�1������������
int point_in_triangle(float(*input1)[4], float(*input2)[4], int m, int n, int k)
{
	float A_x = *(*(input1 + m));
	float A_y = *(*(input1 + m) + 1);
	float B_x = *(*(input1 + n));
	float B_y = *(*(input1 + n) + 1);
	float C_x = *(*(input1 + k));
	float C_y = *(*(input1 + k) + 1);
	float P_x = **input2;
	float P_y = *(*input2 + 1);

	float signOfTrig = (B_x - A_x)*(C_y - A_y) - (B_y - A_y)*(C_x - A_x);

	float signOfAB = (B_x - A_x)*(P_y - A_y) - (B_y - A_y)*(P_x - A_x);
	float signOfCA = (A_x - C_x)*(P_y - C_y) - (A_y - C_y)*(P_x - C_x);
	float signOfBC = (C_x - B_x)*(P_y - C_y) - (C_y - B_y)*(P_x - C_x);

	int r1 = signOfTrig*signOfAB < 0 ? 0 : 1;
	int r2 = signOfTrig*signOfCA < 0 ? 0 : 1;
	int r3 = signOfTrig*signOfBC < 0 ? 0 : 1;

	if ((r1&r2&r3) && signOfTrig*signOfAB*signOfCA*signOfBC)
		return 1;
	else if ((r1&r2&r3) && !(signOfTrig*signOfAB*signOfCA*signOfBC))
		return 0;
	else
		return -1;
}

int point_in_rectangle(float(*input1)[4], float(*input2)[4])
{
	int flag = 0;//flag ����һ���������� |=0x01������һ���������� |=0x02��
	int ret = 0;
	int r1, r2, r3, r4;

	r1 = point_in_triangle(input1, input2, 0, 1, 2);
	r2 = point_in_triangle(input1, input2, 0, 1, 3);
	r3 = point_in_triangle(input1, input2, 0, 2, 3);
	r4 = point_in_triangle(input1, input2, 1, 2, 3);

	if (r1 == 1)//Trig1
	{
		ret |= 0x01;
		flag |= 0x01;
	}
	else if (r1 == 0)
	{
		if (!(flag & 0x02))//�㲻���κ���������
		{
			ret |= 0x01;
			flag |= 0x02;
		}
		else if (!(flag & 0x01) && (flag & 0x02))//����һ����������
			ret |= 0x01;
	}

	if (r2 == 1)//Trig2
	{
		ret |= 0x02;
		flag |= 0x01;
	}
	else if (r2 == 0)
	{
		if (!(flag & 0x02))//�㲻���κ���������
		{
			ret |= 0x02;
			flag |= 0x02;
		}
		else if (!(flag & 0x01) && (flag & 0x02))
			ret |= 0x02;
	}

	if (r3 == 1)//Trig3
	{
		ret |= 0x04;
		flag |= 0x01;
	}
	else if (r3 == 0)
	{
		if (!(flag & 0x02))//�㲻���κ���������
		{
			ret |= 0x04;
			flag |= 0x02;
		}
		else if (!(flag & 0x01) && (flag & 0x02))
			ret |= 0x04;
	}

	if (r4 == 1)//Trig4
	{
		ret |= 0x08;
		flag |= 0x01;
	}
	else if (r4 == 0)
	{
		if (!(flag & 0x02))//�㲻���κ���������
		{
			ret |= 0x08;
			flag |= 0x02;
		}
		else if (!(flag & 0x01) && (flag & 0x02))
			ret |= 0x08;
	}

	return ret;
}

int return_value(float(*input1)[4], float(*input2)[4], float(*output)[4])
{
	int ret = 0;
	int i = 0, j = 0;
	int flag = 0;

	ret = point_in_rectangle(input1, input2);

	if (ret & 0x01)//0,1,2
	{
		if (!flag)
		{
			for (i = 0; i < 4; i++)
				*(*(output + 0) + i) = *(*(input1 + 0) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 1) + i) = *(*(input1 + 1) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 2) + i) = *(*(input1 + 2) + i);
			flag = 1;
		}
		else
		{
			for (i = 0; i < 4; i++)
				*(*(output + 3) + i) = *(*(input1 + 0) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 4) + i) = *(*(input1 + 1) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 5) + i) = *(*(input1 + 2) + i);
		}
	}

	if (ret & 0x02)//0,1,3
	{
		if (!flag)
		{
			for (i = 0; i < 4; i++)
				*(*(output + 0) + i) = *(*(input1 + 0) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 1) + i) = *(*(input1 + 1) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 2) + i) = *(*(input1 + 3) + i);
			flag = 1;
		}
		else
		{
			for (i = 0; i < 4; i++)
				*(*(output + 3) + i) = *(*(input1 + 0) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 4) + i) = *(*(input1 + 1) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 5) + i) = *(*(input1 + 3) + i);
		}
	}

	if (ret & 0x04)//0,2,3
	{
		if (!flag)
		{
			for (i = 0; i < 4; i++)
				*(*(output + 0) + i) = *(*(input1 + 0) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 1) + i) = *(*(input1 + 2) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 2) + i) = *(*(input1 + 3) + i);
			flag = 1;
		}
		else
		{
			for (i = 0; i < 4; i++)
				*(*(output + 3) + i) = *(*(input1 + 0) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 4) + i) = *(*(input1 + 2) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 5) + i) = *(*(input1 + 3) + i);
		}
	}

	if (ret & 0x08)//1,2,3
	{
		if (!flag)
		{
			for (i = 0; i < 4; i++)
				*(*(output + 0) + i) = *(*(input1 + 1) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 1) + i) = *(*(input1 + 2) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 2) + i) = *(*(input1 + 3) + i);
			flag = 1;
		}
		else
		{
			for (i = 0; i < 4; i++)
				*(*(output + 3) + i) = *(*(input1 + 1) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 4) + i) = *(*(input1 + 2) + i);
			for (i = 0; i < 4; i++)
				*(*(output + 5) + i) = *(*(input1 + 3) + i);
		}
	}

	return 0;
}


