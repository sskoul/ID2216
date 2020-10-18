#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>
#include <graphics.h>
#include <fstream>
#include <string>

using namespace std;

double Frobenious_Difference ( vector <vector <double> > A , vector <vector <double> > B) {

	int d = (int)A.size() , i , j;
	
	double Temp = 0;
	
	for ( i = 0 ; i < d; i++ ) {
	
		for ( j = 0; j < d; j++ ) {
		
			Temp += (A[i][j] - B[i][j]) * (A[i][j] - B[i][j]);
		}
	}
	return Temp;
}


//Projecting to Simplex
vector <double> Simplex_Projection (vector <double>& Vec) {
	
	int d = (int)Vec.size(), rho;
	vector <double> z(d) , y(d) , Sum(d) , Res(d);
		
		
	Sum[ 0 ] = Vec[ 0 ];	
	for (int i=1; i < d; i++ ) Sum[ i ] = Sum[ i - 1] + Vec[ i ];
	
	
	for (int i= 0; i < d; i++ ) {
	
			double Temp = double(1) / double(i + 1);
			
			z[i] = Vec[i] + Temp * (1 - Sum[ i ]);\
	}
	

	rho = -1;
	for (int i=0; i < d; i++ ) {
	
			if (z[ i ] > 0 ) rho = i;
			
	}
	
	double lambda = z[rho] - Vec[rho];


	for (int i = 0; i < d; i++) {
	
		Res[i] = max( Vec[i] + lambda , double(0) ) ;
		
	}
	
	return Res;
}

//Projecting to Row Stochastic Matrices
vector <vector <double> > Stochastic_Projection ( vector <vector <double> >& A) {
		
		int d = (int)A.size() , i;
		
		vector <vector <double> > B;
		
		for (i = 0; i < d; i++ ) {
		
			vector <double> Temp = Simplex_Projection( A[ i ] );
			
			B.push_back( Temp );
		}
		
		return B;
}


vector <vector <double> > Transpose( vector <vector <double> >& A ) {

	int d = (int)A.size() , i , j;
	
	vector <vector <double> > B;
	
	vector <double> Temp( d );
	
	for ( i = 0; i < d; i++ ) B.push_back( Temp );

	for ( i = 0; i < d; i++ ) {
	
		for (j =0; j < d; j++ ) {
		
			B[ i ][ j ] = A[ j ][ i ];
		}
		
	}
	
	return B;
}


//Projecting to Column Stochastic Matrices
vector <vector <double> > Col_Stochastic_Projection ( vector <vector <double> >& A) {
		
		vector <vector <double> > B = Transpose( A );
		
		vector <vector <double> > C = Stochastic_Projection( B );
		
		return Transpose( C );
}


//Projecting to Doubly Stochastic Matrices
vector <vector <double> > Doubly_Stochastic_Projection ( vector <vector <double> >& A) {
		
		vector <vector <double> > B = A , C;
		
		int d = (int)A.size();		
		
		
		for (int i=0; i < 3; i++) {
		
			C = Stochastic_Projection(B);
			
			B = Col_Stochastic_Projection( C );
		}
		
		return B;
}



//Computing the gradients of Min Sum Set Cover Cost for a given request
vector <double> Gradient( vector <int> Request, vector <vector <double> >& A ) {

	int element, d = (int)A.size() , i, i_star , S = (int)Request.size();
	
	vector <double> Sum( d );
		
	Sum[ 0 ] = 0;
	
	for ( i = 1; i < d; i++ ) {
	
		double Temp = 0;
		
		for (int k = 0; k < S; k++ ) {

			Temp += A[ Request[ k ] ][ i - 1];
		}
		
		Sum[ i ] = Sum[ i - 1] + Temp;
	}
	 
	for ( i = 0; i < d; i++ ) {
	
		if ( Sum[ i ] > 1 ) {
		
			i_star = i;
			
			break;
		}
	}
	
	if( Sum[d - 1] <= 1 ) i_star = d;
	
	vector <double> Grad( d );
	
	for (i = 0; i < i_star; i++ ) Grad[ i ] =  ( i_star - i );
	
	return Grad;
}



//Simulates Projected Gradient Descent
void Projected_Gradient_Descent( vector <int> Request, vector <vector <double> >& A, int T ) {

	int d = (int)A.size() , element , S = (int)Request.size() , i;
	
	vector <double> grad = Gradient( Request, A );
	

	for (int m = 0; m < S; m++ ) {
	
		element = Request[ m ];
		
		for (int i= 0; i < d; i++ ) {
		
			A[element][ i ] = A[element][ i ] + grad[i] / double( d * d * sqrt(T)); 
		}
	}
		
	A = Doubly_Stochastic_Projection( A );
}



//Step 1 of Algorithm 4 in Appendix B.2
double Sample (vector <double> Prob_Distrib) {

	int Num_Buckets = (int)Prob_Distrib.size();
	
	vector <double> Sum(Num_Buckets);
	

	double r = (double)rand() / (double)RAND_MAX;
		
	Sum[ 0 ]  = Prob_Distrib[ 0 ];
	for (int i = 1; i < Num_Buckets; i++ ) {
			Sum[ i ] = Sum[ i - 1] + Prob_Distrib[ i ];
	}
	
	for(int i =0 ; i < Num_Buckets; i++ ) {
		
		if (r < Sum[i] ) return double(i+1)/(double)Num_Buckets;
	
	}
	
	return 1;
}


vector< vector <double> > Mul_Matrix(vector< vector <double> >&  A , double  a ) {
	int m = (int)A.size(), n = (int)A[0].size();		
	
	vector <vector <double> > B;
	
	for (int i = 0; i < m; i++ ) {

		vector <double> Temp(n);
		
		for (int j = 0; j < n; j++) Temp[ j ] = a * A[ i ][ j ];
	
		B.push_back( Temp );
	}
	return B;
}
 
//Step 1 of Algorithm 4 in Appendix B.2
double Sample_Pareto( int Num_Buckets) {

	vector <double> P (Num_Buckets);
	
	for (int i = 0; i < Num_Buckets; i++ ) {
		
		P[ i ] = double( 2 * i + 1) / double(Num_Buckets * Num_Buckets); 
	
	}
	
	return Sample( P );
}


//Algorithm 4 of Appendix B.2
vector <int> Rand_Rounding( vector <vector <double> >& A ) {
	vector <int> Res;
	int d = (int)A.size();

	double a = Sample_Pareto(100) , Q = 1.67/a;
	
   
    vector <vector <double> > B  = Mul_Matrix( A , Q );
	vector <pair< int,int > > Pairs;
	
	for (int e = 0; e < d; e++ ) {
		for( int j = 0; j <= d/2 - 1; j++ ) {
			B[e][2 * j + 1] = B[ e ][2 * j + 1] + B[e][j];
		}
	}
		
		
	for (int e = 0; e < d; e++ ) {
		
		double r = (double)rand() / (double)RAND_MAX;
		
		vector <double> Sum( d );
		
		Sum[ 0 ] = 0;
		for ( int i = 1 ; i < d; i++ ) Sum[ i ] = Sum[ i - 1 ] + B[ e ][ i - 1];
		
		for (int i = d - 1; i >= 0; i--) {
		
			if ( Sum[ i ] < r) {
				
				Pairs.push_back( make_pair ( i , e ) );
				
				break;
			}	
		}
	}
		
	sort( Pairs.begin() , Pairs.end() );
	
	pair <int,int> p;
	
	for ( int i = 0; i < (int)Pairs.size(); i++ ) {
		
		p = Pairs[ i ];
		
		Res.push_back( p .second);
	}
	
	return Res;	
}

//generate requests
vector<vector <int> > Generate_Requests(int Set_Size, int Num_Elements, int Requests_Num) {

	vector <vector <int> > Requests;
	
	for (int i = 1; i <= Requests_Num; i++ ) {
		
		vector <int> Temp;
		
	
	// uncomment the lines 326-331 and comment lines 333-338 to switch between the types of requests described in Section 5.
	
	//	Temp.push_back( rand() % 2);

	//	for (int j = 1; j <= Set_Size - 1; j++) {
		
	//		Temp.push_back( rand() % Num_Elements) ;
	//	}
	
		Temp.push_back( rand() % 5);

		for (int j = 1; j <= Set_Size - 1; j++) {
		
			Temp.push_back( rand() % Num_Elements) ;
		}
	
	
		
		Requests.push_back(Temp);
	}
	
	return Requests;
}
 
 
//outputs the cost of given request in an given permutation
int Cost ( vector <int> Per, vector <int> Request ) {

	int d = (int)Per.size() , S = (int)Request.size();
	
	for (int i = 0; i < d; i++ ) {
	
		for ( int j = 0; j < S; j++ ) {
			
			if ( Request[ j ] == Per[ i ] ) return i + 1;
		}
	}
	return d;
}

vector <int> Simulate_Rand_Round( vector <vector <int> >& Reqs, int Num_Elements) {

	int d = Num_Elements;

	vector <vector < double>  > A;
	vector <double> Temp;
	vector <int> Per , Res;
	
	for (int i = 1; i <= d; i++ ) Temp.push_back( (double)1/(double)d );
	
	for (int i = 1; i <= d; i++ ) A.push_back( Temp );
	
	for (int l = 0; l < (int)Reqs.size(); l++) {
		
		Per = Rand_Rounding( A );	
		Res.push_back( Cost( Per , Reqs[ l ] ));	
		Projected_Gradient_Descent( Reqs[ l ], A , l + 1 );

	}

	return Res;
}

//Producing a permutation generated uniformly at random
vector <int> Generate_Random_Permutation( int Num_Elements) {

	int d = Num_Elements;
	vector <pair <double , int> > Pairs;
	vector <int> Res;

	
	for (int i = 0; i < d; i++ ) {
	
		double r = (double)rand() / (double)RAND_MAX;
			
		Pairs.push_back( make_pair( r , i ) );
		
	}
	
	sort( Pairs.begin() , Pairs.end() );
	
	pair <double , int> p;
	
	for (int i = 0; i < d; i++ ) {
		
		p = Pairs[ i ];
		
		Res.push_back( p.second );
	}
	
	return Res;

}


//Simulating the cost of randomly selected permutation
vector <int>  Simulate_Rand_Permutation( vector <vector <int> >& Reqs, int Num_Elements) {

	int d = Num_Elements;
	
	vector <int> Per , Res;
	
	
	for (int l = 0; l < (int)Reqs.size(); l++) {
				
		Per = Generate_Random_Permutation( Num_Elements );
		Res.push_back( Cost( Per , Reqs[ l ] ) );
	}
	 
	return Res;
}


void Remove_Element( vector <int>& Vec , int Target) {

	vector <int> Res;
	
	for (int i = 0; i < (int)Vec.size(); i++ ) {
	
		if( Vec[ i ] != Target ) Res.push_back( Vec[ i ] );
	}
	
	Vec = Res;
}


bool hasElement( vector <int> Request , int Target ) {

	for (int i = 0; i < (int)Request.size(); i++ ) {
		
		if ( Request[ i ] == Target ) return true;
	}
	
	return false;
}

// Algorithm of Feige-Lovasz-Tetali
int Find_Most_Frequent( vector <vector <int> >& Requests , int Num_Elements, set <int>& Previously_Placed ) {

	int d  = Num_Elements;
	vector <int> Freqs (d);
	vector <int> Temp;
	set <int> :: iterator it;
	
	for (int i = 0; i < (int)Requests.size(); i++ ) {
	
		Temp = Requests[ i ];
		
		for (int j = 0; j < (int)Temp.size(); j++ ) {
		
			Freqs[ Temp[ j ] ]++;
		}
	}
	
	int Res = -1;
		
	for (int i = 0; i < d; i++ ) {
		
		it = Previously_Placed.find( i );
			
		if (it == Previously_Placed.end() ) {
			
			Res = i;	
			break;
		}
	}
			
	for (int i = 0; i < d; i++ ) {
		
		it = Previously_Placed.find( i );
		
		if (it == Previously_Placed.end() && Freqs[ i ] > Freqs[ Res ] ) Res = i;
	}
		
	return Res;
}

// Algorithm of Feige-Lovasz-Tetali
vector <int> LTF ( vector <vector <int> >& Reqs , int Num_Elements ) {

	vector <int> Res;
	int Temp , d = Num_Elements;
	set <int> Placed;
	
	vector <vector <int> > Requests = Reqs;
	vector < vector <int> > Reqs_Temp;
	
	for (int i = 0; i < d; i++ ) {
	
		Temp = Find_Most_Frequent( Requests , Num_Elements , Placed);
			
		Res.push_back( Temp ); Placed.insert( Temp );
		
		for (int j = 0; j < (int)Requests.size(); j++ ) {
		
			if ( !hasElement( Requests[ j ] , Temp )  )  Reqs_Temp.push_back( Requests[ j ] );
		}
		
		Requests = Reqs_Temp;
		
		Reqs_Temp.erase( Reqs_Temp.begin() ,  Reqs_Temp.end() );
	}
	
	return Res;
}


// Algorithm of Feige-Lovasz-Tetali
vector <int> Simulate_LTF ( vector <vector <int> >& Reqs , int Num_Elements ) {

	vector <int> Per = LTF( Reqs , Num_Elements );
	 
	vector <int> Res;
	
	for (int i = 0; i < (int)Reqs.size(); i++ ) {
		
		Res.push_back( Cost(Per , Reqs[ i ] ) ); 
	}
	
	return Res;
}



void Write_File(string Name , int Number , vector <int>& Vec) {

	char c = '0' + Number;
	
	Name.push_back(c); Name.append(".txt");
	
	char* path = const_cast<char*>(Name.c_str());
    std::ofstream file( path );
    
    for (int i = 0; i < Vec.size(); i++ ) {
	
		file << Vec[i] << "\n";
	}
}


void Sub( vector <double> &A , vector <double>& B ) {

	for (int i = 0; i < (int)A.size(); i++ ) A[ i ] = max( A[ i ] - B[ i ] , (double)0 );

}

// Algorithm 5 of Appendix B.6
int First (vector <double> Target_Vector , vector<vector <double> >& Vecs , set <int>& Prohibited ) {

	int N = (int) Target_Vector.size() , i , j , Res;
	bool flag = false;
	double Value , Min_Val;
	
	if ( (int)Prohibited.size() == N ) return -1;
	
	set <int> :: iterator it;
	
	vector <double> Values( N );
	
	for ( i = 0; i < N; i++ ) {
		
		it = Prohibited.find( i );
		
		if (it != Prohibited.end() ) continue;
	
		Value = 0;

		for ( j = 0; j < N; j++ )  Value += max ( Target_Vector[ j ] - Vecs[ i ][ j ] , (double)0 );

		if ( !flag ) {
		
			Min_Val = Value;
			
			flag = true;
			
			Res = i;
		
		} else {
					if ( Value < Min_Val) {
					
						Min_Val = Value;
						
						Res = i;
					}
			}		
	}
	Prohibited.insert( Res );
	return Res;
}

// Algorithm 5 of Appendix B.6
vector <int> First_r (vector<vector <double> >& Vecs , set <int>& Prohibited, int r ) {

	vector <int> Res;
	vector <double> Target_Vector;
	int N = (int)Vecs[0].size() , Temp;
	
	
	for (int i = 0; i < N; i++ ) Target_Vector.push_back( (double)1 );

	for (int i = 1; i <= r; i++ ) {
	
		Temp = First( Target_Vector , Vecs , Prohibited );
		
		if (Temp == -1) break;
		
		Res.push_back( Temp );
		Sub( Target_Vector , Vecs[ Temp ] );
	}
	
	return Res;
}

void Append( vector <int>& A , vector <int>& B ) {

	for ( int i = 0; i < (int)B.size(); i++ ) A.push_back( B[ i ] );

}

// Algorithm 2 of Section 4
vector <int> Det_Rounding ( vector <vector <double> >& A , int r ) {

	int N = (int)A.size();
	set <int> Prohibited;
	vector <vector <double> > B( N );
	vector <int> Temp (N) , Res;
	vector <double> Temp2( N );
	
	
	for (int i = 0; i < N; i++ ) {
	
		Temp2[ 0 ] = (double)0;
		
		for (int j = 1; j < N; j++ ) Temp2[ j ] = A[ i ][ j - 1] + Temp2[ j - 1 ]; 
		
		B[ i ] = Temp2;
	}
		
	for ( int l = 1; l <= N / r + 1; l++ ) {

		Temp = First_r( B , Prohibited , r);
				
		Append( Res , Temp );	
	}

	return Res;
}


// Simulation of Online Projected Gradient Descent + Deterministic Rounding
vector <int> Simulate_Det_Round( vector <vector <int> >& Reqs, int Num_Elements , int Cardinality) {

	int d = Num_Elements;

	vector <vector < double>  > A;
	vector <double> Temp;
	vector <int> Per , Res;
	
	for (int i = 1; i <= d; i++ ) Temp.push_back( (double)1/(double)d );
	
	for (int i = 1; i <= d; i++ ) A.push_back( Temp );
	
	for (int l = 0; l < (int)Reqs.size(); l++) {
		
		Per = Det_Rounding( A , Cardinality);	
		Res.push_back( Cost( Per , Reqs[ l ] ));	
		Projected_Gradient_Descent( Reqs[ l ], A , l + 1 );

	}

	return Res;
}



// Simulation of Deterministic , Randomized , Feige-Lovasz-Tetali. The argument "Times" denote the number of times the simulation is conducted and its value must agree
//with the value of the variable "Copies" in "Make_Plots.py".

void Repeat (int Rounds, int Cardinality, int Num_Elements, int Times) {

	vector <int> Cost_Rand_Round , Cost_Rand_Per , Cost_LTF , Cost_Det_Round;

	for (int i = 1; i <=Times; i++ ) {
	
		vector <vector <int> > Reqs = Generate_Requests( Cardinality , Num_Elements , Rounds );
			
			Cost_Rand_Round = Simulate_Rand_Round( Reqs , Num_Elements );

			Cost_Det_Round = Simulate_Det_Round( Reqs , Num_Elements , Cardinality );
	
			Cost_Rand_Per = Simulate_Rand_Permutation( Reqs , Num_Elements );
	
			Cost_LTF = Simulate_LTF( Reqs , Num_Elements );
	
			Write_File( "Random_Permutation" , i , Cost_Rand_Per  );
			Write_File( "Det_Rounding" , i , Cost_Det_Round  );
			Write_File( "Random_Rounding" , i , Cost_Rand_Round  );
			Write_File( "LTF" , i , Cost_LTF  );
	}
	
	return;
}



int main () {
	srand (time(0));

	//Repeat(10000 , 5 , 100 , 3 );

	Repeat(10000 , 10 , 100 , 3 );
	return 1;
}

