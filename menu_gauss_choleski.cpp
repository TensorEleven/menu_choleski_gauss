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
    void loadMat(string filename);
    void decompMat();            //retourne une matrice triangulaire sup B
    void transpose(float **mat); //transposer une matrice
	void solveTriangSup(float** mat, float* sec);       //resolution de matrice triangulaire superieure
    void solveTriangInf();       //resolution de matrice triangulaire inferieure
    
    void process(int choice);

    //method de resolution
    void solveCholeski();
    void solveGauss();
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
    //float *resTriangSup;
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

/// remplissage des données
    
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
    fichier.close();
}

void Lsolver::loadMat(string filename){
    size_t i(0), j(0);
    ifstream fichier(filename, ios::in);
    if(fichier){
        fichier >> dim;

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
    fichier.close();
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
            float sum = 0;           //somme des B[i][k].B[j][k]
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

//Résoudre une matrice triangulaire supérieur tq
// mat.res = y
// mat : matrice carré d'ordre dim
// res : 
void Lsolver::solveTriangSup(float** mat,float* sec){ //mat=Bt, sec=y
    float s(0);
    int i(0), j(0);
    for(i=dim-1; i>=0; i--){    /// Must go backward
        for(j=i+1, s=0; j<int(dim); j++)
            s += (mat[i][j]*x[j]);
        x[i] = (sec[i]-s)/mat[i][i];
    }
}

void Lsolver::solveTriangInf(){ //mat=B, sec=b
    float s(0);
    int i(0), j(0);
    for(i=0; i<(int)dim; i++){  /// Must go forward
        for(j=0, s=0; j<int(dim); j++)
            s += (B[i][j]*y[j]);
        y[i] = (b[i]-s)/B[i][i];
    }
}
/*
void Lsolver::solveTriangInf(){ //mat=B, sec=b, res=y
    float s(0);
    int i(0), j(0);
    for(i=0; i<(int)dim; i++){  /// Must go forward
        for(j=0, s=0; j<int(dim); j++)
            s += (B[i][j]*y[j]);
        y[i] = (b[i]-s)/B[i][i];
    }
}*/

void Lsolver::solveCholeski(){
    /*
     *  Choleski
     */
	decompMat();
    //display
    cout << "\nDécomposistion de la matrice A tq A = B.Bt" << endl;
    cout << "\nLa matrice B :" << endl;
    displayMat(dim,B);

    cout << "\nLa matrice Bt :" << endl;
    displayMat(dim,Bt);
    
    //définir les valeurs de y en tand que solution du systeme B.y = b
    solveTriangInf();
    //solver.solveTriangInf();

    //définir les valeurs de x en tant que solution du system Bt.x = y
    solveTriangSup(Bt,y);
}

void Lsolver::solveGauss(){//m=A, btmp=b
    //repeter l'operation jusqu'à dim-1
    for (int k=0; k<(int)dim;k++){
        //trouver pivot max
        float piv = abs(A[k][k]);
        int lpiv = k;
        
        for(int i=k;i<(int)dim;i++){
            if(piv<abs((int)A[i][k])){
                piv = abs(A[i][k]);   //trouver le max |pivot|
                lpiv = i;                   //sauvegarder la igne
            }
        }
        //echanger la ligne k+1 au ligne lpiv
        for (int i=0; i<(int)dim; i++){
            std::swap<float>(A[k][i],A[lpiv][i]);
        }
        std::swap<float>(b[k],b[lpiv]);
       
        //reduire par method de pivot
        for(int i=k+1;i<(int)dim;i++){
            for (int j=k+1; j<(int)dim; j++){  //or for (int j=k;j<dim;j++)
                if(A[i][j]==0)
                    continue;
                A[i][j] = A[i][j] - (A[i][k]/A[k][k])*A[k][j];
            }
            b[i] = b[i] - (A[i][k]/A[k][k])*b[k];
            A[i][k]=0;
        }
        //trouver la solution
        solveTriangSup(A,b);
    }
}

/*
 * 
 *  M E N U
 * 
 * 
 */
void Lsolver::process(int choice) {
    switch (choice)
    {
    case 1:
        cout << "Method de Gauss" << endl;
        solveGauss();
        break;
        
    case 2:
        cout << "Choleski" << endl;
        solveCholeski();
        break;
    default:
        break;
    }
}
// Structure de données pour le menu
struct Menu {
    Menu() {
        menuItemsLength = 2;
        menuItems = nullptr;
        menuItems = new (nothrow) string[menuItemsLength];
        menuItems[0] = "Résoudre par méthode de Gauss";
        menuItems[1] = "Résoudre par méthode de Choleski";
    }

    int display() {
        int i = 0;
        cout << "Choisi un operation :" << endl;
        for(i = 0; i < menuItemsLength; i++) {
            cout << "\t" << i + 1 << " : " << menuItems[i] << endl;
        }
        cout << "\t0" << " : Quit" << endl;
        cout << "\t> : "; cin >> choice;
        return choice;
    }

    string * menuItems;
    int menuItemsLength;
    int choice;
};

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
    
/// Resolution
    //solver.solveCholeski();
    //solver.solveGauss();
    Menu menu;
    int choice=10;
    while(true) {
        choice = menu.display();
        if (choice == 0) {
            break;
        }
        solver.process(choice);

        /// Résultats
	    solver.displayResult();
        solver.loadMat("data.txt");
    }

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
