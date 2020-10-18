import matplotlib.pyplot as plt
import numpy as np
import math


def convert ( s ):

	Res = 0
	
	for i in range( 0 , len(s) - 1) :

		Bit = ord( s[ i ]) - ord('0')
		
		Res = 10 * Res + Bit

	return Res


def Time_Average ( A ):

	Res = [];
	
	Res.append ( A[ 0 ] );
	
	for i in range( 1 , len(A) ):

		Res.append( Res[ i - 1 ] + A[ i ])

	Res2 = []

	for i in range( 0 , len(A) ):

		Res2.append( Res[ i ] / (i + 1) )
		
	return Res2


def expected_vector ( A ):

	Temp = A[ 0 ]
	
	Res = np.asarray( Temp )
	
	for i in range( 1 , len(A) ):
		
		Res = Res + np.asarray( A[ i ] )
	
	return Res / len( A )


def variance_vector ( A ):

	Vector_Len = len( A[ 0 ] )

	Exp = expected_vector( A )
	
	Res = []

	for i in range( 0 , Vector_Len ):
	
		Temp = 0
	
		for j in range(0 , len( A )):
		
			Temp = Temp + ( A[ j ][ i ] - Exp[i] )* ( A[ j ][ i ] - Exp[i] )

		Res.append( math.sqrt( (Temp / len( A ) )))

	return np.asarray( Res )
	
	
def Open_Files( Name , Num_Copies ):

	Res = []

	for l in range(1 , Num_Copies + 1):

		c = chr( 48 + l)
	
		f = open(Name + c + ".txt", 'r');
		data = f.readlines()
		
		Temp = []
	
		for i in range(0 , len(data) ):
	
			Temp.append( convert( data[ i ] ) )
	
		Res.append(Time_Average( Temp) ) ;

	return Res

Copies = 3

Random_Permutation = Open_Files( "Random_Permutation" , Copies )
Det_Rounding = Open_Files( "Det_Rounding" , Copies )
Random_Rounding = Open_Files( "Random_Rounding" , Copies )
LTF = Open_Files( "LTF" , Copies )

EXP_Random_Permutations = expected_vector( Random_Permutation )
VAR_Random_Permutations = variance_vector( Random_Permutation )

EXP_Det_Rounding = expected_vector( Det_Rounding )
VAR_Det_Rounding = variance_vector( Det_Rounding )

EXP_Random_Rounding = expected_vector( Random_Rounding )
VAR_Random_Rounding = variance_vector( Random_Rounding )

EXP_LTF = expected_vector( LTF )
VAR_LTF = variance_vector( LTF )

Len = len(EXP_LTF )

plt.figure(figsize=(10,8));

plt.ylabel('Total Access Cost / #Requests')
plt.xlabel('#Requests')	

plt.plot(EXP_Random_Permutations, color='red', linewidth=3, label=r'Random Permutation');
plt.plot(EXP_Random_Permutations + VAR_Random_Permutations,color='tab:red');
plt.plot(EXP_Random_Permutations - VAR_Random_Permutations,color='tab:red');
plt.fill_between(np.arange(0, Len),EXP_Random_Permutations - VAR_Random_Permutations,EXP_Random_Permutations + VAR_Random_Permutations,color='tab:red', alpha=1)

plt.plot(EXP_Det_Rounding, color='purple', linewidth=3 , label=r'OPDG with Deterministic Rounding');
plt.plot(EXP_Det_Rounding + VAR_Det_Rounding,color='tab:purple');
plt.plot(EXP_Det_Rounding - VAR_Det_Rounding,color='tab:purple');
plt.fill_between(np.arange(0, Len),EXP_Det_Rounding - VAR_Det_Rounding,EXP_Det_Rounding + VAR_Det_Rounding,color='tab:purple', alpha=1)

plt.plot(EXP_Random_Rounding, color='blue', linewidth=3, label=r'OPDG with Randomized Rounding');
plt.plot(EXP_Random_Rounding + VAR_Random_Rounding,color='tab:blue');
plt.plot(EXP_Random_Rounding - VAR_Random_Rounding,color='tab:blue');
plt.fill_between(np.arange(0, Len),EXP_Random_Rounding - VAR_Random_Rounding,EXP_Random_Rounding + VAR_Random_Rounding,color='tab:blue', alpha=1)


plt.plot(EXP_LTF, color='green', linewidth=3, label=r'4-approx (offline) by Lovasz Feige Tetali');
plt.plot(EXP_LTF, color='green', linewidth=3);
plt.plot(EXP_LTF + VAR_LTF,color='tab:green');
plt.plot(EXP_LTF - VAR_LTF,color='tab:green');
plt.fill_between(np.arange(0, Len),EXP_LTF - VAR_LTF,EXP_LTF + VAR_LTF,color='tab:green', alpha=1)

plt.text(2000,26,'#Elements = 100',horizontalalignment='center')
plt.text(2000,25,'Cardinality of Requests = 5',horizontalalignment='center')

plt.tick_params(labelsize=10)
plt.legend(fontsize=10)
plt.savefig('plot.png', dpi=300, bbox_inches='tight')

#plt.show()