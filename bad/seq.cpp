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
	double pos_x, pos_y, vel_x, vel_y, acc_x, acc_y, ast_m;
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
			ast_m = 0; // Masa
		}
		// Setters
		void set_pos_x(double n){pos_x = n;}
		void set_pos_y(double n){pos_y = n;}
		void set_vel_x(double n){vel_x = n;}
		void set_vel_y(double n){vel_y = n;}
		void set_acc_x(double n){acc_x = n;}
		void set_acc_y(double n){acc_y = n;}
		void set_m(double n){ast_m = n;}
		// Getters
		double get_pos_x(){return pos_x;}
		double get_pos_y(){return pos_y;}
		double get_vel_x(){return vel_x;}
		double get_vel_y(){return vel_y;}
		double get_acc_x(){return acc_x;}
		double get_acc_y(){return acc_y;}
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
 * Funciones:
 *	calcDistAst
 *	calcFuerzaAstX
 *	calcFuerzaAstY
 *	calcDistPla
 *	calcFuerzaPlaX
 *	calcFuerzaPlaY
 *	checkArgs
 */
/**
 * calcDistAst
 * Calcula la distancia entre el Asteroide a y el Asteroide b a partir de sus posiciones.
 * Retorna un double.
 */
double calcDistAst(Asteroide a, Asteroide b)
{
	return sqrt(pow(a.get_pos_x() - b.get_pos_x(), 2) + pow(a.get_pos_y() - b.get_pos_y(), 2));
}
/**
 * calcFuerzaAstX
 * Calcula la fuerza ejercida sobre el Asteroide a por el Asteroide b en el eje X.
 * Hace uso de la constante de gravitación universal, g.
 * Retorna un double
 */
double calcFuerzaAstX(double g, Asteroide a, Asteroide b)
{
	double distancia = calcDistAst(a, b);
	double pendiente = (a.get_pos_y() - b.get_pos_y()) / (a.get_pos_x() - b.get_pos_x());
	if (pendiente > 1 || pendiente < -1)
	{
		pendiente -= trunc(pendiente);
	}
	double angulo = atan(pendiente);
	return g * a.get_m() * b.get_m() * cos(angulo) / pow(distancia, 2);
}
/**
 * calcFuerzaAstY
 * Calcula la fuerza ejercida sobre el Asteroide a por el Asteroide b en el eje Y.
 * Hace uso de la constante de gravitación universal, g.
 * Retorna un double
 */
double calcFuerzaAstY(double g, Asteroide a, Asteroide b)
{
	double distancia = calcDistAst(a, b);
	double pendiente = (a.get_pos_y() - b.get_pos_y()) / (a.get_pos_x() - b.get_pos_x());
	if (pendiente > 1 || pendiente < -1)
	{
		pendiente -= trunc(pendiente);
	}
	double angulo = atan(pendiente);
	return g * a.get_m() * b.get_m() * sin(angulo) / pow(distancia, 2);
}
/**
 * calcDistPla
 * Calcula la distancia entre el Asteroide a y el Planeta b a partir de sus posiciones.
 * Retorna un double.
 */
double calcDistPla(Asteroide a, Planeta b)
{
	return sqrt(pow(a.get_pos_x() - b.get_pos_x(), 2) + pow(a.get_pos_y() - b.get_pos_y(), 2));
}
/**
 * calcFuerzaPlaX
 * Calcula la fuerza ejercida sobre el Asteroide a por el Planeta b en el eje X.
 * Hace uso de la constante de gravitación universal, g.
 * Retorna un double
 */
double calcFuerzaPlaX(double g, Asteroide a, Planeta b)
{
	double distancia = calcDistPla(a, b);
	double pendiente = (a.get_pos_y() - b.get_pos_y()) / (a.get_pos_x() - b.get_pos_x());
	if (pendiente > 1 || pendiente < -1)
	{
		pendiente -= trunc(pendiente);
	}
	double angulo = atan(pendiente);
	return g * a.get_m() * b.get_m() * cos(angulo) / pow(distancia, 2);
}
/**
 * calcFuerzaPlaY
 * Calcula la fuerza ejercida sobre el Asteroide a por el Planeta b en el eje Y.
 * Hace uso de la constante de gravitación universal, g.
 * Retorna un double
 */
double calcFuerzaPlaY(double g, Asteroide a, Planeta b)
{
	double distancia = calcDistPla(a, b);
	double pendiente = (a.get_pos_y() - b.get_pos_y()) / (a.get_pos_x() - b.get_pos_x());
	if (pendiente > 1 || pendiente < -1)
	{
		pendiente -= trunc(pendiente);
	}
	double angulo = atan(pendiente);
	return g * a.get_m() * b.get_m() * sin(angulo) / pow(distancia, 2);
}
/**
 * checkArgs
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
	initFile = fopen("init_conf.txt", "w");
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
	for (int i = 0; i < num_iteraciones; ++i)
	{
		// printf("Iteración %d\n", i);
		for (unsigned int j = 0; j < vec_ast.size(); ++j) // Iterador principal de Asteroides
		{
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
			double fuerza_x = 0, fuerza_y = 0; // Valores acumulativos de la fuerza en cada eje
			//Cálculo de fuerzas Asteroide -> Asteroide
			for (unsigned int k = 0; k < vec_ast.size(); ++k) // Iterador secundario de Asteroides
			{
				if (k != j) // No hay efecto sobre sí mismo
				{
					if (calcDistAst(vec_ast[j], vec_ast[k]) > dmin)
					{
						// Acumulación de fuerzas en las variables fuerza_x y fuerza_y
						fuerza_x += calcFuerzaAstX(gravity, vec_ast[j], vec_ast[k]);
						fuerza_y += calcFuerzaAstY(gravity, vec_ast[j], vec_ast[k]);
					}
				}
			}
			// Cálculo de fuerzas Planeta -> Asteroide
			for (unsigned int k = 0; k < vec_pla.size(); ++k) // Iterador de Planetas
			{
				// Acumulación de fuerzas en las variables fuerza_x y fuerza_y
				fuerza_x += calcFuerzaPlaX(gravity, vec_ast[j], vec_pla[k]);
				fuerza_y += calcFuerzaPlaY(gravity, vec_ast[j], vec_pla[k]);
			}
			// Almacenamiento de la aceleración en cada eje para cada asteroide
			// !!! No se aplican los efectos hasta calcular las aceleraciones de todos los asteroides !!!
			vec_ast[j].set_acc_x(fuerza_x / vec_ast[j].get_m());
			vec_ast[j].set_acc_y(fuerza_y / vec_ast[j].get_m());
		}
		// Aplicación de las aceleraciones para obtener las nuevas posiciones y velocidades
		for (unsigned int j = 0; j < vec_ast.size(); ++j)
		{
			vec_ast[j].set_vel_x(vec_ast[j].get_vel_x() + vec_ast[j].get_acc_x() * delta_t);
			vec_ast[j].set_vel_y(vec_ast[j].get_vel_y() + vec_ast[j].get_acc_y() * delta_t);
			vec_ast[j].set_pos_x(vec_ast[j].get_pos_x() + vec_ast[j].get_vel_x() * delta_t);
			vec_ast[j].set_pos_y(vec_ast[j].get_pos_y() + vec_ast[j].get_vel_y() * delta_t);
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
	outFile = fopen("out.txt", "w");
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
