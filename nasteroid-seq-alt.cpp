#include <iostream>
#include <random>
#include <stdlib.h>
#include <math.h>
using namespace std;

/**
 * Clases:
 *	Asteroide
 *	Planeta
 *	Rayo
 */
class Asteroide
{
	double pos_x, pos_y, vel_x, vel_y, acc_x, acc_y, fue_x, fue_y, ast_m;
	public:
		// Constructor vacío
		Asteroide()
		{
			pos_x = 0; // Coordenada x
			pos_y = 0; // Coordenada y
			vel_x = 0; // Velocidad x
			vel_y = 0; // Velocidad y
			acc_x = 0; // Aceleración x
			acc_y = 0; // Aceleración y
			fue_x = 0; // Fuerza x
			fue_y = 0; // Fuerza y
			ast_m = 0; // Masa
		}
		// Setters
		void set_pos_x(double n){pos_x = n;}
		void set_pos_y(double n){pos_y = n;}
		void set_vel_x(double n){vel_x = n;}
		void set_vel_y(double n){vel_y = n;}
		void set_acc_x(double n){acc_x = n;}
		void set_acc_y(double n){acc_y = n;}
		void set_fue_x(double n){fue_x = n;}
		void set_fue_y(double n){fue_y = n;}
		void set_m(double n){ast_m = n;}
		// Getters
		double get_pos_x(){return pos_x;}
		double get_pos_y(){return pos_y;}
		double get_vel_x(){return vel_x;}
		double get_vel_y(){return vel_y;}
		double get_acc_x(){return acc_x;}
		double get_acc_y(){return acc_y;}
		double get_fue_x(){return fue_x;}
		double get_fue_y(){return fue_y;}
		double get_m(){return ast_m;}
};

class Planeta
{
	double pos_x, pos_y, pla_m;
	public:
		// Constructor vacío
		Planeta()
		{
			pos_x = 0; // Coordenada x
			pos_y = 0; // Coordenada y
			pla_m = 0; // Masa
		}
		// Setters
		void set_pos_x(double n){pos_x = n;}
		void set_pos_y(double n){pos_y = n;}
		void set_m(double n){pla_m = n;}
		// Getters
		double get_pos_x(){return pos_x;}
		double get_pos_y(){return pos_y;}
		double get_m(){return pla_m;}
};
class Rayo
{
	double pos_x, pos_y;
	public:
		// Constructor no vacío
		Rayo(double x, double y)
		{
			pos_x = x; // Coordenada x (0)
			pos_y = y; // Coordenada y
		}
		// Setters
		void set_pos_x(double n){pos_x = n;}
		void set_pos_y(double n){pos_y = n;}
		// Getters
		double get_pos_x(){return pos_x;}
		double get_pos_y(){return pos_y;}		
};

/**
 * Función: checkArgs
 * Comprueba la validez de los argumentos de entrada del programa.
 * Debe haber 6 argumentos en total (1 + 5).
 * Deben ser convertibles a int con la función stoi.
 * Retorna 0 en caso positivo o detiene el programa en caso negativo.
 */
int checkArgs(int argc, char *argv[])
{
	if (argc != 6)
	{
		cout << "\nnasteroids-seq: Wrong arguments.\n";
		cout << "Correct use:\n";
		cout << "nasteroids-seq num_asteroides num_iteraciones num_planetas pos_rayo semilla\n\n";
		exit (-1);
	}
	for (int i = 1; i < argc; ++i)
	{
		try
		{
			stoi(argv[i]);
		}
		catch (invalid_argument& e)
		{
			cout << "\nnasteroids-seq: Wrong arguments.\n";
			cout << "Correct use:\n";
			cout << "nasteroids-seq num_asteroides num_iteraciones num_planetas pos_rayo semilla\n\n";
			exit(-1);
		}
	}
	return 0;
}

/**
 * main
 */
int main(int argc, char *argv[])
{
	checkArgs(argc, argv); // Se validan los argumentos antes de nada

	// Se almacenan los argumentos en variables con nombres más cómodos
	int num_asteroides	= stoi(argv[1]);
	int num_iteraciones	= stoi(argv[2]);
	int num_planetas	= stoi(argv[3]);
	double pos_rayo		= stod(argv[4]);
	int semilla			= stoi(argv[5]);

	// Se definen las constantes del problema
	const double gravity	= 6.674e-5;
	const double delta_t	= 0.1;
	const double dmin		= 2.0;
	const double width		= 200;
	const double height		= 200;
	const double ray_width	= 4;
	const double mass		= 1000;
	const double sdm		= 50;

	// Función generadora de distribuciones aleatorias
	default_random_engine re{semilla};
	uniform_real_distribution<double>xdist{0.0,std::nextafter(width, std::numeric_limits<double>::max())};
	uniform_real_distribution<double>ydist{0.0,std::nextafter(height, std::numeric_limits<double>::max())};
	normal_distribution<double>mdist{mass,sdm};

	// Vector que almacena instancias de la clase Asteroide
	vector <Asteroide> vec_ast;
	for (int i = 0; i < num_asteroides; ++i)
	{
		vec_ast.push_back(Asteroide());		// Creación del objeto Asteroide
		vec_ast[i].set_pos_x(xdist(re));	// Asignación de la coordenada x
		vec_ast[i].set_pos_y(ydist(re));	// Asignación de la coordenada y
		vec_ast[i].set_m(mdist(re));		// Asignación de la masa
	}

	// Vector que almacena instancias de la clase Planeta
	vector <Planeta> vec_pla;
	for (int i = 0; i < num_planetas; ++i)
	{
		switch (i%4) // Hay 4 formas de crear los Planetas. Se distinguen mediante la operación módulo 4.
		{
			case 0:	vec_pla.push_back(Planeta());		// Creación del objeto Planeta
					vec_pla[i].set_pos_x(0);			// Asignación de la coordenada x
					vec_pla[i].set_pos_y(ydist(re));	// Asignación de la coordenada y
					vec_pla[i].set_m(mdist(re) * 10);	// Asignación de la masa
					break;
			case 1:	vec_pla.push_back(Planeta());
					vec_pla[i].set_pos_x(xdist(re));
					vec_pla[i].set_pos_y(0);
					vec_pla[i].set_m(mdist(re) * 10);
					break;
			case 2:	vec_pla.push_back(Planeta());
					vec_pla[i].set_pos_x(width);
					vec_pla[i].set_pos_y(ydist(re));
					vec_pla[i].set_m(mdist(re) * 10);
					break;
			case 3:	vec_pla.push_back(Planeta());
					vec_pla[i].set_pos_x(xdist(re));
					vec_pla[i].set_pos_y(height);
					vec_pla[i].set_m(mdist(re) * 10);
					break;
		}
	}

	Rayo death_star(0, pos_rayo); // Creación del objeto Rayo y asignación de las coordenadas

	// Se almacenan los datos generados en el fichero init_conf.txt
	FILE * initFile;
	initFile = fopen("init_conf-alt.txt", "w");
	fprintf(initFile, "%d %d %d %.3f %d\n", num_asteroides, num_iteraciones, num_planetas, pos_rayo, semilla); // Argumentos introducidos
	for (unsigned int i = 0; i < vec_ast.size(); ++i)
	{
		fprintf(initFile, "%.3f %.3f %.3f\n", vec_ast[i].get_pos_x(), vec_ast[i].get_pos_y(), vec_ast[i].get_m()); // Vector de Asteroides
	}
	for (unsigned int i = 0; i < vec_pla.size(); ++i)
	{
		fprintf(initFile, "%.3f %.3f %.3f\n", vec_pla[i].get_pos_x(), vec_pla[i].get_pos_y(), vec_pla[i].get_m()); // Vector de Planetas
	}
	fprintf(initFile, "%.3f %.3f\n", death_star.get_pos_x(), death_star.get_pos_y()); // Rayo
	fclose(initFile);

	// Bucle para el número de iteraciones
	// Variables temporales utilizadas en los cálculos, sólo se declaran una vez.
	double distancia_no_sqrt, distancia_real;
	double pendiente, angulo, temp_fuer, temp_fuer_x, temp_fuer_y, ast1_pos_x, ast1_pos_y, ast1_mas, ast2_pos_x, ast2_pos_y, ast2_mas, pla1_pos_x, pla1_pos_y, pla1_mas;
	for (int i = 0; i < num_iteraciones; ++i)
	{
		for (unsigned int j = 0; j < vec_ast.size(); ++j) // Iterador principal de Asteroides
		{
			// Se extraen los datos del asteroide j
			ast1_pos_x	= vec_ast[j].get_pos_x();	// Asteroide j, coordenada x
			ast1_pos_y	= vec_ast[j].get_pos_y();	// Asteroide j, coordenada y
			ast1_mas	= vec_ast[j].get_m();		// Asteroide j, masa

			//Cálculo de fuerzas Asteroide -> Asteroide
			for (unsigned int k = 0; k < vec_ast.size(); ++k) // Iterador secundario de Asteroides
			{
				// Se extraen los datos del asteroide k
				ast2_pos_x	= vec_ast[k].get_pos_x();	// Asteroide k, coordenada x
				ast2_pos_y	= vec_ast[k].get_pos_y();	// Asteroide k, coordenada y
				ast2_mas	= vec_ast[k].get_m();		// Asteroide k, masa

				// Distancia sin aplicar la raíz cuadrada, para más precisión en el cálculo de la fuerza, ya que se evita hacer la raíz y después elevar al cuadrado
				distancia_no_sqrt = ((ast1_pos_x - ast2_pos_x) * (ast1_pos_x - ast2_pos_x)) + ((ast1_pos_y - ast2_pos_y) * (ast1_pos_y - ast2_pos_y));
				distancia_real = sqrt(distancia_no_sqrt);

				if (distancia_real > dmin)
				{
					pendiente = (ast2_pos_y - ast1_pos_y) / (ast2_pos_x - ast1_pos_x); // Cálculo de la pendiente
					if (!isinf(pendiente) && (pendiente > 1 || pendiente < -1))
					{
						pendiente -= trunc(pendiente); // Truncado de la pendiente
					}
					angulo = atan(pendiente); // Cálculo del ángulo de efecto

					temp_fuer = gravity * ast1_mas * ast2_mas / distancia_no_sqrt; // Fuerza total
					if (temp_fuer > 200)
					{
						temp_fuer = 200; // Truncado de la fuerza
					}

					temp_fuer_x = temp_fuer * cos(angulo); // Componente x de la fuerza total
					temp_fuer_y = temp_fuer * sin(angulo); // Componente y de la fuerza total

					// vec_ast[j].set_fue_x(vec_ast[j].get_fue_x() - temp_fuer_x); // Al asteroide j se le suma
					// vec_ast[j].set_fue_y(vec_ast[j].get_fue_y() - temp_fuer_y);
					vec_ast[k].set_fue_x(vec_ast[k].get_fue_x() + temp_fuer_x); // Al asteroide k se le resta
					vec_ast[k].set_fue_y(vec_ast[k].get_fue_y() + temp_fuer_y);
				}
			}
		}
		for (unsigned int j = 0; j < vec_pla.size(); ++j)
		{
			pla1_pos_x	= vec_pla[j].get_pos_x();
			pla1_pos_y	= vec_pla[j].get_pos_y();
			pla1_mas	= vec_pla[j].get_m();

			// Cálculo de fuerzas Planeta -> Asteroide
			for (unsigned int k = 0; k < vec_ast.size(); ++k) // Iterador de Planetas
			{
				// Acumulación de fuerzas en las variables fuerza_x y fuerza_y
				ast2_pos_x	= vec_ast[k].get_pos_x();
				ast2_pos_y	= vec_ast[k].get_pos_y();
				ast2_mas	= vec_ast[k].get_m();

				distancia_no_sqrt = ((pla1_pos_x - ast2_pos_x) * (pla1_pos_x - ast2_pos_x)) + ((pla1_pos_y - ast2_pos_y) * (pla1_pos_y - ast2_pos_y));

				pendiente = (ast2_pos_y - pla1_pos_y) / (ast2_pos_x - pla1_pos_x);
				if (!isinf(pendiente) && (pendiente > 1 || pendiente < -1))
				{
					pendiente -= trunc(pendiente);
				}
				angulo = atan(pendiente);

				temp_fuer = gravity * pla1_mas * ast2_mas / distancia_no_sqrt;
				if (temp_fuer > 200)
				{
					temp_fuer = 200;
				}

				temp_fuer_x = temp_fuer * cos(angulo);
				temp_fuer_y = temp_fuer * sin(angulo);

				vec_ast[k].set_fue_x(vec_ast[k].get_fue_x() + temp_fuer_x);
				vec_ast[k].set_fue_y(vec_ast[k].get_fue_y() + temp_fuer_y);
			}
		}
		// Aplicación de las aceleraciones para obtener las nuevas posiciones y velocidades
		for (unsigned int j = 0; j < vec_ast.size(); ++j)
		{
			// Se calculan y aplican las aceleraciones
			vec_ast[j].set_acc_x(vec_ast[j].get_fue_x() / vec_ast[j].get_m());
			vec_ast[j].set_acc_y(vec_ast[j].get_fue_y() / vec_ast[j].get_m());

			// Se calculan y aplican las velocidades
			vec_ast[j].set_vel_x(vec_ast[j].get_vel_x() + vec_ast[j].get_acc_x() * delta_t);
			vec_ast[j].set_vel_y(vec_ast[j].get_vel_y() + vec_ast[j].get_acc_y() * delta_t);

			// Se calculan y aplican las nuevas posiciones
			vec_ast[j].set_pos_x(vec_ast[j].get_pos_x() + vec_ast[j].get_vel_x() * delta_t);
			vec_ast[j].set_pos_y(vec_ast[j].get_pos_y() + vec_ast[j].get_vel_y() * delta_t);

			// Se ponen las fuerzas a 0 para la siguiente iteración.
			vec_ast[j].set_fue_x(0);
			vec_ast[j].set_fue_y(0);

			// El efecto rebote se aplica después de los cálculos, para asegurar que no se guardan valores erróneos en out.txt
			if (vec_ast[j].get_pos_x() <= 0) // Efecto rebote x <= 0
			{
				vec_ast[j].set_pos_x(2);
				vec_ast[j].set_vel_x(-vec_ast[j].get_vel_x());
			}
			if (vec_ast[j].get_pos_y() <= 0) // Efecto rebote y <= 0
			{
				vec_ast[j].set_pos_y(2);
				vec_ast[j].set_vel_y(-vec_ast[j].get_vel_y());
			}
			if (vec_ast[j].get_pos_x() >= width) // Efecto rebote x >= width
			{
				vec_ast[j].set_pos_x(width - 2);
				vec_ast[j].set_vel_x(-vec_ast[j].get_vel_x());
			}
			if (vec_ast[j].get_pos_y() >= height) // Efecto rebote y >= height
			{
				vec_ast[j].set_pos_y(height - 2);
				vec_ast[j].set_vel_y(-vec_ast[j].get_vel_y());
			}
		}
		// Ejecución del rayo destructor
		// Bucle que recorre el vector de asteroides y elimina los que deban ser destruidos con erase()
		// Credit: https://stackoverflow.com/a/23383451
		if(!vec_ast.empty())
		{
			for(int j = vec_ast.size() - 1; j >= 0; j--)
			{
				if(abs(vec_ast.at(j).get_pos_y() - death_star.get_pos_y()) <= ray_width/2)
				{
					vec_ast.erase(vec_ast.begin() + j); 
				}
			}
		}
	}

	// Se almacenan los resultados en el fichero out.txt
	FILE * outFile;
	outFile = fopen("out-alt.txt", "w");
	for (unsigned int i = 0; i < vec_ast.size(); ++i)
	{
		fprintf(outFile, "%.3f %.3f %.3f %.3f %.3f\n",	vec_ast[i].get_pos_x(),
														vec_ast[i].get_pos_y(),
														vec_ast[i].get_vel_x(),
														vec_ast[i].get_vel_y(),
														vec_ast[i].get_m());
	}
	fclose(outFile);

	return 0;
}
