/*
 * R�solution de syst�me lin�aire par M�thode choleski
 * V 1.0
 * 2021-06-22
 * Auteur(s): ROBINSON Olivier, RATEFIARISON Harivony
 */

#include <iostream>
#include <fstream>
#include <cmath>
//#include <locale>
#include <iomanip>

using namespace std;

/// General helpers
float **newMat(int rows, int cols);
template <class T>
	T *newVect(size_t dim);
void displayMat(size_t dim, float **A);
void displayVec(size_t dim, float *v);

class Lsolver{
public:
	Lsolver(string filename);    //filename est nom du fichier contenant les données de la matrice
	~Lsolver();
	void displayResult();
    void decompMat();            //retourne une matrice triangulaire sup B
    void transpose(float **mat); //transposer une matrice
	void solveTriangSup();       //resolution de matrice triangulaire superieure
    void solveTriangInf();       //resolution de matrice triangulaire inferieure
    
	size_t  getdim(){ return dim;}
	float** getMat(){ return A;}
    float** getB(){return B;}
    float** getBt(){return Bt;}
	float*  getRhs(){ return b;}
private:
	size_t dim;
	float  **A;		    // matrice du probl�me A.x=b
    float **B;          // matrice triangulaire sup d udecomposition
	float **Bt;         // matrice transposé de B
    float  *b, *x, *y;	// second membre et inconnu du probl�me
};

Lsolver::Lsolver(string filename){
/// ouverture du fichier de donn�es
	size_t i(0), j(0);
    ifstream fichier(filename, ios::in);
    if(fichier){
        fichier >> dim;
/// allocation de la matrice, du second membre, B, Bt, y et de la solution
        A = newMat(dim, dim);
        b = newVect<float> (dim);   //vecteur formé du second membre
        B = newMat(dim,dim);        //matrice de decomposition, triang sup
        Bt = newMat(dim,dim);       //matrice transposee de B
        x = newVect<float> (dim);   //solution du systeme
        y = newVect<float> (dim);   //vecteur tq B.y = b

/// remplissage des donn�es
        for(i=0; i<dim; i++){
            for(j=0; j<dim; j++){
                fichier >> A[i][j]; 
            }
        }
        for(i=0; i<dim; i++){
            fichier >> b[i];
		}
        for(i=0; i<dim; i++){       //initialiser les valeurs des vecteurs
            x[i] = 0;
            y[i] = 0;
		}
        fichier.close();            //fermer le fichier
    }
    else{
        cout << "Données non trouvées..." << endl;  //message d'erreur
    }
}

//dectructeur pour les tableaux alloués

Lsolver::~Lsolver(){
	delete[] x;
	delete[] y;
    delete[] b;

	for(size_t i=0; i< dim; i++){
        delete[] A[i];
        delete[] B[i];
        delete[] Bt[i];
    }
	delete[] A;
    delete[] B;
    delete[] Bt;
}

void Lsolver::displayResult(){
	//float eps(1e-6);
    cout << "\nVecteur y:" << endl;
    for(size_t i=0; i<(size_t)dim; i++)
        cout << "y" << i+1 << " = " << y[i] << endl;
    cout << "\nLa solution" << endl;
    for(size_t i=0; i<dim; i++)
        cout << "x" << i+1 << " = " << x[i] << endl;
}


void Lsolver::decompMat(){           //trouver B,Bt tq A=B.Bt et A.x=b
    for(int i=0; i<(int)dim;i++){
        for(int j=0; j<=i;j++){
            float sum = 0;          //somme des B[i][k].B[j][k]
            if(i!=j){
                for(int k=0;k<=j-1;k++){
                    sum+=(B[i][k]*B[j][k]);
                }
                B[i][j] = (1/B[j][j]*(A[i][j]-sum));
            }
            else if (i==j){       //i==j
                for(int k=0;k<=i-1;k++){
                    sum += pow(B[i][k],2);
                }
                B[i][i]= sqrt(A[i][i]- sum);
            }
        }
    }
    transpose(B);                       //remplir Bt avec le transposé de B
}

void Lsolver::transpose(float **mat){   //remplir Bt avec le transposé de la matrice mat
    for(int i=0;i<(int)dim;i++){
        for(int j=0;j<(int)dim;j++){
            Bt[i][j] = mat[j][i];
        }
    }
}

void Lsolver::solveTriangSup(){
    float s(0);
    int i(0), j(0);
    for(i=dim-1; i>=0; i--){    /// Must go backward
        for(j=i+1, s=0; j<int(dim); j++)
            s += (Bt[i][j]*x[j]);
        x[i] = (y[i]-s)/Bt[i][i];
    }
}
void Lsolver::solveTriangInf(){
    float s(0);
    int i(0), j(0);
    for(i=0; i<(int)dim; i++){  /// Must go forward
        for(j=0, s=0; j<int(dim); j++)
            s += (B[i][j]*y[j]);
        y[i] = (b[i]-s)/B[i][i];
    }
}

int main(){
    //Afficher l'environnement
	//cout << "Locale: " << setlocale(LC_ALL,"") << endl;
    cout << "Résolution d'un système linéaire A.x=b" << endl;
    
/// Données
    Lsolver solver("data.txt");
    cout << "La matrice A : " << endl;
    displayMat(solver.getdim(), solver.getMat());
    cout << "Le second membre b : " << endl;
    displayVec(solver.getdim(), solver.getRhs());
    
/// Calculs
	solver.decompMat();     //
    //display
    cout << "\nDécomposistion de la matrice A tq A = B.Bt" << endl;
    cout << "\nLa matrice B :" << endl;
    displayMat(solver.getdim(),solver.getB());

    cout << "\nLa matrice Bt :" << endl;
    displayMat(solver.getdim(),solver.getBt());
    solver.solveTriangInf();
    solver.solveTriangSup();
/// R�sultats
	solver.displayResult();

    return 0;
}

void displayMat(size_t dim, float **A){
	cout << "[";
	for(size_t i=0; i<dim; i++){
		if(i==0)cout << "[";
		else    cout << " [";
		for(size_t j=0; j<dim; j++){
			if(j<dim-1)cout << setw(8) << setprecision(5) << A[i][j] << "  ";
			else cout << setw(8) << setprecision(5) << A[i][j];
		}
		if(i<dim-1) cout << "]" << endl;
		else 		cout << "]]" << endl;
	}
}

void displayVec(size_t dim, float *v){
	for(size_t i=0; i<dim; i++){
		cout << " [" << v[i] << "]" << endl;
	}
}

template <class T>
T *newVect(size_t dim){
    T *v(NULL);
    v = new T[dim];
    if(v==NULL){
        cout << "Erreur lors de la création du vecteur" << endl;
    }
    return v;
}

float **newMat(int rows, int cols){
    float **mat(NULL);
    mat = newVect<float *>(rows);

    for(int i=0; i<rows; i++)
        mat[i] = newVect<float>(cols);
    return mat;
}
