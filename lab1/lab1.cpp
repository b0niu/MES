#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

#define DOKL 0.00000000001

class GlobalData
{

public:

	double H;
	double W;
	double nH;
	double nW;
	double k;
	int np;
	double ro;
	double c;
	double alfa;

	double Simt;
	double dt;

	double Tot;
	double Tpocz;
	
	double nN;
	double nE;
	double dy;
	double dx;

	GlobalData() {
		double* dane = Dane(13);

		H = dane[0];
		W = dane[1];
		nH = dane[2];
		nW = dane[3];
		k = dane[4];
		np = (int)dane[5];
		c = dane[6];
		ro = dane[7];
		alfa = dane[8];

		Simt = dane[9];
		dt = dane[10];

		Tot = dane[11];
		Tpocz = dane[12];

		nN = dane[3] * dane[2];
		nE = (dane[3] - 1) * (dane[2] - 1);
		dy = dane[0] / (dane[2] - 1);
		dx = dane[1] / (dane[3] - 1);
	}

	double* Dane(int it) {
		double* d = new double[it];
		fstream dane;
		dane.open("dane.txt", ios::in);

		for (int i = 0; i < it; i++)
		{
			dane >> d[i];
		}

		dane.close();

		return d;
	}

};

class NODE {
public:
	int ind;
	double x, y;
	double T;
	bool BC;
	NODE(int i, double X, double Y, bool bc, double t) {
		ind = i;
		x = X;
		y = Y;
		T = t;
		BC = bc;
	}
	NODE() {
		x = 0.0;
		y = 0.0;
		T = 0.0;
		ind = 1;
		BC = false;
	}
public:
	friend ostream& operator<< (ostream& wyjscie, const NODE& n)
	{
		return wyjscie << fixed << setprecision(3) << "(x->" << n.x << " y->" << n.y << ")" << "[" << n.T << "]";
		
	}
	void* operator new(size_t rozmiar) {
		void* wsk = malloc(rozmiar);
		return wsk;
	}

	
};

class Uklad_Row
{
	double** A;
	double* X;
	double* F;
	int n;
public:
	Uklad_Row(double** a, double* f, double* x,int size) {
		A = a;
		F = f;
		X = x;
		n = size;

		if (rozklad_LU(n, A) && wynik_LU(n, A, F, X))
		{
			for (int i = 0; i < n; i++)
			{
				//cout << "x" << i + 1 << " = " << X[i] << endl;
			}
		}
		else  cout << "Dzielenie przez zero!!!" << endl;
		

	}
	~Uklad_Row() {

	}

	bool rozklad_LU(int n, double** A)
	{
		for (int i = 0; i < n - 1; i++)
		{
			//Gdy dzielimy przez zero -> Koniec
			if (fabs(A[i][i]) < DOKL) return false;

			for (int j = i + 1; j < n; j++)
			{
				A[j][i] /= A[i][i];
			}

			for (int j = i + 1; j < n; j++)
			{
				for (int k = i + 1; k < n; k++)
				{
					A[j][k] -= A[j][i] * A[i][k];
				}
			}
		}
		return true;
	}

	bool wynik_LU(int n, double** A, double* F, double* X)
	{
		double suma; //Suma iloczynów

		X[0] = F[0];

		for (int i = 1; i < n; i++)
		{
			suma = 0;
			for (int j = 0; j < i; j++)
			{
				suma += A[i][j] * X[j]; //Wyliczamy i aktualizujemy sume
			}
			X[i] = F[i] - suma;
		}

		//Gdy dzielimy przez zero -> Koniec
		if (fabs(A[n - 1][n - 1]) < DOKL) return false;

		X[n - 1] /= A[n - 1][n - 1]; // Wyznaczamy i obliczamy wektor X

		for (int i = n - 2; i >= 0; i--)
		{
			suma = 0;
			for (int j = i + 1; j < n; j++)
			{
				suma += A[i][j] * X[j];
			}
			//Gdy dzielimy przez zero -> Koniec
			if (fabs(A[i][i]) < DOKL) return false;

			X[i] = (X[i] - suma) / A[i][i]; // Wyznaczamy wartości X
		}
		return true;
	}

private:

};



class ELEM {
public:
	NODE id[4];
	int ind;

	ELEM(int I) {
		ind = I;
	}

	ELEM(){
		ind = 1;
	}

	int Count() {
		return 4;
	}

	int countID(int i) {

		return i;
	}


	void setN(NODE N1, NODE N2, NODE N3, NODE N4) {
		id[0] = N1;
		id[1] = N2;
		id[2] = N3;
		id[3] = N4;
	}

	NODE* nodes() {
		return id;
	}
	friend ostream& operator<< (ostream &wyjscie, const ELEM& el)
	{
		return wyjscie << el.id[3] << " " << el.id[2] << endl << el.id[0] << " " << el.id[1] << endl;
	}
	void* operator new(size_t rozmiar) {
		void* wsk = malloc(rozmiar);
		return wsk;
	}
	void* operator new[](size_t rozmiar) {
		void* wsk = malloc(rozmiar);
		return wsk;
	}
};

class Grid
{
public:
	GlobalData glData;
	NODE* node;
	ELEM* elements;
	int el_size;

	Grid();
	~Grid();
	void GridsetNodes() {
		node = new NODE[glData.nN];
		int index = 1;
		bool bc;
		for (int i = 0; i < glData.nW; i++) {
			for (int j = 0; j < glData.nH; j++) {
				bc = false;
				if (i * glData.dx == 0 || i * glData.dx == glData.W || j * glData.dy == 0 || j * glData.dy == glData.H) bc = true;
				NODE node1(index, i * glData.dx, j * glData.dy,bc,glData.Tpocz);
				node[index - 1] = node1;
				index++;
			}
		}
	}

	void GridsetNodes_T(double* X) {
		node = new NODE[glData.nN];
		int index = 1;
		bool bc;
		for (int i = 0; i < glData.nW; i++) {
			for (int j = 0; j < glData.nH; j++) {
				bc = false;
				if (i * glData.dx == 0 || i * glData.dx == glData.W || j * glData.dy == 0 || j * glData.dy == glData.H) bc = true;
				NODE node1(index, i * glData.dx, j * glData.dy, bc, X[index-1]);
				node[index - 1] = node1;
				index++;
			}
		}
	}

	void GridsetElem() {

		elements = new ELEM[glData.nE];
		int iE = 1;
		int t1;
		for (int i = 0; i < glData.nE; i++) {
			elements[i] = ELEM(i + 1);
			if ((iE % (int)glData.nH == 0) && (i != 0)) {
				iE++;
			}
			t1 = node[iE - 1].ind - 1;
			elements[i].setN(
				node[t1],
				node[node[iE - 1].ind - 1 + (int)glData.nH],
				node[node[iE - 1].ind - 1 + (int)glData.nH + 1],
				node[node[iE - 1].ind - 1 + 1]);
			iE++;
			
		}

	}

	void showMax_Min() {
		double min = node[1].T;
		double max = node[1].T;
			for (int i = 0; i < glData.nN; i++) {
				if (min > node[i].T)min = node[i].T;
				if (max < node[i].T)max = node[i].T;
			}
			cout << fixed;
			cout << setprecision(3) << "\tMax T= " << max << endl << "\tMin T= " << min << endl;
			fstream plik;
			fstream plik1;
			plik.open("max2.txt", ios::app);
			plik1.open("min2.txt", ios::app);
			plik << max << endl;
			plik1 << min << endl;

			plik.close();
			plik1.close();

	}
	ELEM* getElem() { return elements; }

	void changeT(double* X) {
		
		GridsetNodes_T(X);
		GridsetElem(); 
		showMax_Min();
	}

};

Grid::Grid()
{
	GridsetNodes();
	GridsetElem();

}

Grid::~Grid()
{}


class ElLocal
{
public:

	double* pc;
	double* w = new double[2];
	double* NdKsi;
	double* NdEta;
	double* N;
	ElLocal(double ksi, double eta, double w1, double w);
	~ElLocal();


};

ElLocal::ElLocal(double ksi, double eta, double w1, double w2) {
	pc = new double[2]{ ksi,eta };
	w[0] = w1;
	w[1] = w2;
	double N1 = 0.25 * (1 - ksi) * (1 - eta);
	double N2 = 0.25 * (1 + ksi) * (1 - eta);
	double N3 = 0.25 * (1 + ksi) * (1 + eta);
	double N4 = 0.25 * (1 - ksi) * (1 + eta);

	double N1Ksi = (-0.25) * (1 - eta);
	double N2Ksi = 0.25 * (1 - eta);
	double N3Ksi = 0.25 * (1 + eta);
	double N4Ksi = (-0.25) * (1 + eta);

	double N1Eta = (-0.25) * (1 - ksi);
	double N2Eta = (-0.25) * (1 + ksi);
	double N3Eta = 0.25 * (1 + ksi);
	double N4Eta = 0.25 * (1 - ksi);

	N = new double[4]{N1, N2, N3, N4};
	NdKsi = new double[4] {N1Ksi, N2Ksi, N3Ksi, N4Ksi};
	NdEta = new double[4] {N1Eta, N2Eta, N3Eta, N4Eta};
}

ElLocal::~ElLocal()
{
}

class MATRIX {
public:
	//współrzędne pkt całk
	double pc1 = (1.0 / sqrt(3));
	double w1 = 1;

	double pc2_1 = (sqrt(3.0 / 5.0));
	double pc2_2 = 0.0;
	double w2_1 = (5.0 / 9.0);
	double w2_2 = (8.0 / 9.0);

	double pc3_1 = (sqrt((3.0 / 7.0) - 2.0 / 7.0 * sqrt(6.0 / 5.0)));
	double pc3_2 = (sqrt((3.0 / 7.0) + 2.0 / 7.0 * sqrt(6.0 / 5.0)));
	double w3_1 = ((18.0 + sqrt(30)) / 36.0);
	double w3_2 = ((18.0 - sqrt(30)) / 36.0);


	GlobalData gldata;
	 double K = gldata.k;
	int PC = gldata.np;
	double c = gldata.c;
	double ro = gldata.ro;
	ElLocal* local;
	
	double* KSI;
	double* W;

	double** N;
	double**  NdKsi;
	double** NdEta;
	double** J;
	double* detJ;
	double** Ndy;
	double** Ndx;
	
	double** Bc;
	double* Pl;
	double** H;
	double** C;
	double** glC;
	double** glH;
	double* glP;

	int llength;

	MATRIX() {
		switch (PC) {
		case 2: {
			local = new ElLocal[4]{
					ElLocal(-pc1, -pc1, w1, w1),
					ElLocal(pc1, -pc1, w1, w1),

					ElLocal(pc1, pc1, w1, w1),
					ElLocal(-pc1, pc1, w1, w1)
			};
			KSI = new double[2]{ -pc1,pc1 };
			W = new double[2]{ w1,w1 };
			llength = 4;
		}
			  break;
		case 3: {
			local = new ElLocal[9]{
					ElLocal(-pc2_1, -pc2_1, w2_1, w2_1),
					ElLocal(pc2_2, -pc2_1, w2_2, w2_1), 
					ElLocal(pc2_1, -pc2_1, w2_1, w2_1),

					ElLocal(-pc2_1, pc2_2, w2_1, w2_2),
					ElLocal(pc2_2, pc2_2, w2_2, w2_2),
					ElLocal(pc2_1, pc2_2, w2_1, w2_2),

					ElLocal(-pc2_1, pc2_1, w2_1, w2_1),
					ElLocal(pc2_2, pc2_1, w2_2, w2_1),
					ElLocal(pc2_1, pc2_1, w2_1, w2_1)
			};
			KSI = new double[3]{ -pc2_1 ,pc2_2, pc2_1 };
			W = new double[3]{ w2_1,w2_2, w2_1 };
			llength = 9;
		}
			  break;
		case 4: {
			local = new ElLocal[16]{
					ElLocal(-pc3_2, -pc3_2, w3_2, w3_2),
					ElLocal(-pc3_1, -pc3_2, w3_1, w3_2),
					ElLocal(pc3_1, -pc3_2, w3_1, w3_2),
					ElLocal(pc3_2, -pc3_2, w3_2, w3_2),

					ElLocal(-pc3_2, -pc3_1, w3_2, w3_1),
					ElLocal(-pc3_1, -pc3_1, w3_1, w3_1),
					ElLocal(pc3_1, -pc3_1, w3_1, w3_1),
					ElLocal(pc3_2, -pc3_1, w3_2, w3_1),

					ElLocal(-pc3_2, pc3_1, w3_2, w3_1),
					ElLocal(-pc3_1, pc3_1, w3_1, w3_1),
					ElLocal(pc3_1, pc3_1, w3_1, w3_1),
					ElLocal(pc3_2, pc3_1, w3_2, w3_1),

					ElLocal(-pc3_2, pc3_2, w3_2, w3_2),
					ElLocal(-pc3_1, pc3_2, w3_1, w3_2),
					ElLocal(pc3_1, pc3_2, w3_1, w3_2),
					ElLocal(pc3_2, pc3_2, w3_2, w3_2)
			};
			KSI = new double[4]{ -pc3_2,-pc3_1,pc3_1, pc3_2 };
			W = new double[4]{ w3_2,w3_1, w3_1, w3_2 };
			llength = 16;
		}
			  break;
		}

		
		N = new double*[llength];
		NdKsi = new double* [llength];	
		NdEta = new double*[llength];
		J = new double*[llength];
		detJ = new double[llength];
		Ndy = new double*[llength];
		Ndx = new double*[llength];
		H = new double*[llength];

		for (int i = 0; i < llength; i++) {
			N[i] = new double[4];
			NdKsi[i] = new double[4];
			NdEta[i] = new double[4];
			J[i] = new double [4];
			Ndy[i] = new double[4];
			Ndx[i] = new double[4];
			H[i] = new double[llength];
		
		}
	}

	void matN() {
		for (int i = 0; i < llength; i++) {
			N[i] = new double [4] {
					local[i].N[0],
					local[i].N[1],
					local[i].N[2],
					local[i].N[3]
			};
			
		}
	}

	void matNdKsi() {
		for (int i = 0; i < llength; i++) {
			NdKsi[i] = new double [4] {
					local[i].NdKsi[0],
					local[i].NdKsi[1],
					local[i].NdKsi[2],
					local[i].NdKsi[3]
			};
		}

	}

	void matNdEta() {
		for (int i = 0; i < llength; i++) {
			NdEta[i] = new double [4] {
					local[i].NdEta[0],
					local[i].NdEta[1],
					local[i].NdEta[2],
					local[i].NdEta[3]
			};
		}
	}

	void matJ(ELEM el) {
		double* Xdksi = new double[llength] {0};
		double* Xdeta = new double[llength] {0};
		double* Ydksi = new double[llength] {0};
		double* Ydeta = new double[llength] {0};
		
	

		for (int i = 0; i < llength; i++) {
			for (int j = 0; j < 4; j++) {
				Xdksi[i] += NdKsi[i][j] * el.id[j].x;
				Xdeta[i] += NdEta[i][j] * el.id[j].x;
				Ydksi[i] += NdKsi[i][j] * el.id[j].y;
				Ydeta[i] += NdEta[i][j] * el.id[j].y;
			}
			J[i][0] = Xdksi[i];
			J[i][1] = Xdeta[i];
			J[i][2] = Ydksi[i];
			J[i][3] = Ydeta[i];
			

		}
		
	}

	void det() {

		for (int j = 0; j < llength; j++) {
			detJ[j] = (J[j][0] * J[j][3]) - (J[j][2] * J[j][3]);
			
		}
		
	}

	void setNdx() {
		for (int i = 0; i < llength; i++) {
			for (int j = 0; j < 4; j++) {
				Ndx[i][j] = (1 / detJ[i]) * ((J[i][3] * NdKsi[i][j]) - (J[i][2] * NdEta[i][j]));
				
			}
		}
	}

	void setNdy() {
		for (int i = 0; i < llength; i++) {
			for (int j = 0; j < 4; j++) {
				Ndy[i][j] = (1 / detJ[i]) * ((J[i][3] * NdEta[i][j]) - (J[i][2] * NdKsi[i][j]));
			}
		}
	}

	void matH_C() {

		double*** elMatH = new double** [llength];
		double*** elMatC = new double** [llength];
		for (int i = 0; i < llength; i++) {
			elMatH[i] = new double* [4]{ 0 };
			elMatH[i][0] = new double[4]{ 0 };
			elMatH[i][1] = new double[4]{ 0 };
			elMatH[i][2] = new double[4]{ 0 };
			elMatH[i][3] = new double[4]{ 0 };
			elMatC[i] = new double* [4]{ 0 };
			elMatC[i][0] = new double[4]{ 0 };
			elMatC[i][1] = new double[4]{ 0 };
			elMatC[i][2] = new double[4]{ 0 };
			elMatC[i][3] = new double[4]{ 0 };
		}
		
		H = new double*[4];
		C = new double* [4];
		for (int i = 0; i < 4; i++) {
			H[i] = new double[4]{ 0 };
			C[i] = new double[4]{ 0 };
		}

		for (int i = 0; i < llength; i++) {
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					elMatH[i][j][k] = K * (Ndx[i][j] * Ndx[i][k] + Ndy[i][j] * Ndy[i][k]) * detJ[i] * local[i].w[0] * local[i].w[1];
					elMatC[i][j][k] = c * ro * (N[i][j] * N[i][k] ) * detJ[i] * local[i].w[0]  * local[i].w[1];

					H[j][k] += elMatH[i][j][k];
					C[j][k] += elMatC[i][j][k];
				}
			}
		}
	
	}

	void matBC(Grid grid, int elNumb) {
		Bc = new double* [4];
		Pl = new double[4]{ 0 };
	
		for (int i = 0; i < 4; i++) {
			Bc[i] = new double[4]{ 0 };	
		}
		
		for (int i = 0; i < 4; i++) {

			NODE N1 = grid.elements[elNumb].id[i];
			NODE N2 = grid.elements[elNumb].id[(i + 1) % 4];

			if (!(N1.BC == true && N2.BC == true))
				continue;
			for (int j = 0; j < sqrt(llength); j++) {
				double ksi = 0, eta = 0;
				if (i == 0) {
					ksi = KSI[j];
					eta = -1.0;
				}
				else if (i == 1) {
					eta = KSI[j];
					ksi = 1.0;
				}
				else if (i == 2) {
					ksi = KSI[j];
					eta = 1.0;
				}
				else if (i == 3) {
					eta = KSI[j];
					ksi = -1.0;
				}

				double* N2D = new double[4];

				N2D[0] = 0.25 * (1.0 - ksi) * (1.0 - eta);
				N2D[1] = 0.25 * (1.0 + ksi) * (1.0 - eta);
				N2D[2] = 0.25 * (1.0 + ksi) * (1.0 + eta);
				N2D[3] = 0.25 * (1.0 - ksi) * (1.0 + eta);

				double DETJ = (sqrt(pow(N1.x - N2.x, 2) + pow(N1.y - N2.y, 2)))/2;

				
				for (int k = 0; k < 4; k++) {
					for (int z = 0; z < 4; z++) {
						Bc[k][z] += (N2D[k] * N2D[z] * W[j] * DETJ * gldata.alfa);
						
					}
					
				}
				
				for (int k = 0; k < 4; k++) {
					Pl[k] += (N2D[k] * DETJ *W[j]* gldata.alfa* gldata.Tot);
				
				}
				
			}

		}
		
	}

	void calGlH(Grid grid, double***Hzw) {

		ELEM* glHelem = grid.elements;
		glH = new double* [gldata.nN]{ 0 };
		for (int i = 0; i < gldata.nN; i++) {
			glH[i] = new double[gldata.nN]{ 0 };
		}

		
		for (int i = 0; i < gldata.nE; i++) {
			int* ID = new int[4];
			for (int j = 0; j < 4; j++) {
				ID[j] = glHelem[i].id[j].ind - 1;
			}

			
				for (int z = 0; z < 4; z++) {
					for (int k = 0; k < 4; k++) {
						glH[ID[z]][ID[k]] += Hzw[i][z][k];
					}
				}
			
		}

		/*cout << endl;
		cout << "Matrix H global" << endl;
		for (int i = 0; i < gldata.nN; i++) {
			for (int j = 0; j < gldata.nN; j++) {
				cout << fixed;
				cout << setprecision(3) << glH[i][j] << "\t";
			}
			cout << endl;
		}*/
	}

	void calGlP(Grid grid, double **P) {
		ELEM* glPelem = grid.elements;
		glP = new double[gldata.nN]{ 0 };
		

		for (int i = 0; i < gldata.nE; i++) {
			int* ID = new int[4];
			for (int j = 0; j < 4; j++) {
				ID[j] = glPelem[i].id[j].ind - 1;
			}

			for (int z = 0; z < 4; z++) {
					glP[ID[z]] += P[i][z];
			}
		}

		/*cout << endl;
		cout << "Vector P global" << endl;
		for (int i = 0; i < gldata.nN; i++) {
				cout << fixed;
				cout << setprecision(3) << glP[i] << "\t";
		}
		cout << endl;*/
	}

	void calGlC(Grid grid) {
		ELEM* glHelem = grid.elements;
		glC = new double* [gldata.nN]{ 0 };
		for (int i = 0; i < gldata.nN; i++) {
			glC[i] = new double[gldata.nN]{ 0 };
		}

		for (int i = 0; i < gldata.nE; i++) {
			int* ID = new int[4];
			for (int j = 0; j < 4; j++) {
				ID[j] = glHelem[i].id[j].ind - 1;
			}

			for (int z = 0; z < 4; z++) {
				for (int k = 0; k < 4; k++) {
					glC[ID[z]][ID[k]] += C[z][k] / gldata.dt;
				}
			}
		}

		/*cout << endl;
		cout << "Matrix C global" << endl;
		for (int i = 0; i < gldata.nN; i++) {
			for (int j = 0; j < gldata.nN; j++) {
				cout << fixed;
				cout << setprecision(3) << glC[i][j] << "\t";
			}
			cout << endl;
		}*/
	}

	void calH_C() {
		for (int i = 0; i < gldata.nN; i++) {
			for (int j = 0; j < gldata.nN; j++) {
				glH[i][j] += glC[i][j];
			}
		}

		/*cout << endl;
		cout << "Matrix H+C/dt global" << endl;
		for (int i = 0; i < gldata.nN; i++) {
			for (int j = 0; j < gldata.nN; j++) {
				cout << fixed;
				cout << setprecision(3) << glH[i][j] << "\t";
			}
			cout << endl;
		}*/
	}

	void calP_C(Grid grid) {
		double* T0 = new double[gldata.nN]{ 0 };
		for (int i = 0; i < gldata.nN; i++) {
			T0[i] = grid.node[i].T;
			
		}

		double *temp=new double[gldata.nN]{ 0 };
		//cout << "Vector T0 " << endl;
		for (int i = 0; i < gldata.nN; i++) {
			
			for (int j = 0; j < gldata.nN; j++) {
				
				temp[i]+= glC[i][j] * T0[j];
			}
			
			/*cout << fixed;
			cout << setprecision(3) << T0[i] << "\t";*/
		}


		for (int i = 0; i < gldata.nN; i++) {
			glP[i] += temp[i];
			
		}

		/*cout << endl;
		cout << "Vector P+(c/dt)*T0 global" << endl;
		for (int i = 0; i < gldata.nN; i++) {
			cout << fixed;
			cout << setprecision(3) << glP[i] << "\t";
		}
		cout << endl;*/
	}
};

int main() {
	
	GlobalData g;
	//cout << g.H << " " << g.W << " " << g.k << " " << g.nE << " " << g.nN << " " << g.c << " " << g.ro << " " << g.alfa <<  endl;

	Grid grid;
	MATRIX mac;
	mac.matN();
	mac.matNdKsi();
	mac.matNdEta();
	double*** H= new double**[g.nE];
	double** P = new double* [g.nE];

	
	for (int q = 0; q < g.Simt; q += g.dt) {
	
		for (int z = 0; z < mac.gldata.nE; z++) {
			mac.matJ(grid.getElem()[z]);
			mac.det();
			mac.setNdx();
			mac.setNdy();
			mac.matH_C();
			mac.matBC(grid, z);
			
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					mac.H[i][j] += mac.Bc[i][j];
				}
				H[z] = mac.H;
				P[z] = mac.Pl;
			}
		
		}
			cout << "Czas : " << q+ g.dt << endl;
			mac.calGlH(grid, H);
			mac.calGlC(grid);
			mac.calGlP(grid, P);
			mac.calH_C();
			mac.calP_C(grid);

			double* T1 = new double[g.nN]{ 0 };
			Uklad_Row uk_row(mac.glH, mac.glP, T1, g.nN);
			grid.changeT(T1);
			//cout << endl;

		
	}


	return 0;
}

