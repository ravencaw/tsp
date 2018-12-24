/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package tsp;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class TSP_2 {

	/**
	 * Programa principal donde ejecutamos las implementaciones del problema.
	 * Determinamos el tiempo de ejecución de los algoritmos de cada problema,
	 * almacenandolos en un array para su posterior muestra.
	 * 
	 * 1- BackTracking 2- Busqueda Local 3- Divide y Vencerás
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		System.out.println("\nEXPERIMENTO DEL TSP\n\n");
		int iteraciones = 100;
		Random rnd = new Random();
		// Instancias del problema en los archivos tsp
		String[] instancias = { "berlin52.tsp",
				// "a280.tsp","kroA100.tsp", "kroA150.tsp", "kroA200.tsp",
				// "vm1084.tsp", "vm1748.tsp" , "usa13509.tsp"
		};

		double[][] duraciones = new double[instancias.length][3];

		// Modificar luego para todas las instancias (mira instancias.lenght)
		for (int i = 0; i < instancias.length; i++) {
			System.out.println("\n\n--------------------------------------------");
			System.out.println("Instancia " + (i + 1));
			double[][] dist = leerFichero("data\\" + instancias[i]);
			int nciudades = dist.length;
			double[] s1 = new double[nciudades];
			int[] p = new int[nciudades];
			int repeticiones = 500;

			// //Backtraking
			// System.out.println("Ejecutando Bactracking");
			double start = System.currentTimeMillis();
			// for (int j = 0; j < repeticiones; j++) {
			//// Backtracking(dist, nciudades);
			// }
			double stop = System.currentTimeMillis();
			// duraciones[i][0] = (stop - start) / repeticiones;

			// Busqueda local
			System.out.println("Ejecutando "
					// + "Búsqueda Local"
					+ "Algoritmo Voraz");
			start = System.currentTimeMillis();
			for (int j = 0; j < repeticiones; j++) {
				p =
						// TSPLocalSearch(3, nciudades, 1432, dist,
						GreedyAlgorithm(dist, 4)
				// )
				;
			}
			stop = System.currentTimeMillis();
			duraciones[i][1] = (stop - start) / repeticiones;

			// //Divide y Vencerás
			// System.out.println("Ejecutando Divide y Vencerás");
			// start = System.currentTimeMillis();
			// for (int j = 0; j < repeticiones; j++) {
			// solve(dist, 0);
			// }
			//
			// stop = System.currentTimeMillis();
			// duraciones[i][2] = (stop - start) / repeticiones;
			// System.out.println("--------------------------------------------");
		}

		System.out.println("\n\n\n--------------------------------------------");
		System.out.println("Resultados finales");
		for (int i = 0; i < instancias.length; i++) {
			System.out.println("\nInstancia " + (i + 1));
			System.out.println("\tBT: " + duraciones[i][0]); // BACTRACKING
			System.out.println("\tBL: " + duraciones[i][1]); // BUSQUEDA LOCAL
			System.out.println("\tDV: " + duraciones[i][2]); // DIVIDE Y VENCERAS
		}
	}

	/**
	 * Códico del algoritmo de Backtracking
	 * 
	 * @param dist
	 * @param nciudades
	 * @throws Exception
	 */
	// public static void Backtracking(double[][] dist, int nciudades) throws
	// Exception {
	//
	// double [] solucion = creaMatrizD(nciudades);
	// double[] recorrido = creaMatrizD(nciudades);
	//
	// int cont=0;
	// for (int i = 0; i < dist.length; i++) {
	// recorrido[i] = dist[cont][i];
	//
	// }
	//
	// for (int trabajador = 0; trabajador < solucion.length; trabajador++) {
	// solucion[trabajador] = TSP_2.mejorTarea(dist, recorrido, trabajador);
	// }
	//
	//
	//// permute(recorrido, dist, solucion);
	//
	// System.out.print("\nMejor solucion:");
	// for (int i = 0; i < nciudades; i++) {
	// System.out.print(" " + solucion[i]);
	// }
	//
	//// System.out.println(" Distancia --------> " + evaluate(solucion, dist));
	// }
	/**
	 * Comprueba si esa ciudad ha sido escogida
	 * 
	 * @param x
	 * @param camino
	 * @param ciudad
	 * @return
	 * @throws Exception
	 */

	/**
	 * Divide y venceras
	 * 
	 * @param entrada-
	 * @return
	 */
	public static double[][] divideYVenceras(double[][] entrada) {
		double[][][] subS = new double[2][][], subP;
		double[][] s;
		if (entrada.length == 1) {
			s = entrada;
		} else {
			subP = divide(entrada);
			for (int i = 0; i < 2; i++) {
				subS[i] = divideYVenceras(subP[i]);
			}
			s = combina(subS, entrada.length);
		}
		return s;
	}

	/**
	 * Busqueda local en Java Algoritmo que se encargara de buscar la mejor solucion
	 * posible
	 * 
	 * @param vecindad
	 * @param num_ciudades
	 * @param seed
	 * @param dist
	 * @param solucion
	 * @return
	 */
	public static int[] TSPLocalSearch(int vecindad, int num_ciudades, int seed, double[][] dist, int[] solucion) {

		boolean fallo = false;
		int[] auxiliar = creaMatriz(num_ciudades); // para evitar los fallos de punteros al trabajar con la solucion
		System.arraycopy(solucion, 0, auxiliar, 0, solucion.length);

		while (!fallo && vecindad >= 0) { // mientras no haya fallo y la vecindad sea mayor o igual a 0 (para la
											// recursividad)
			auxiliar = cambioCiudades(auxiliar, num_ciudades, seed);
			if (evaluate(auxiliar, dist) < evaluate(solucion, dist)) {
				for (int j = 0; j < num_ciudades; j++) { // Para copiar un array a otro
					solucion[j] = auxiliar[j];
				}
				System.arraycopy(TSPLocalSearch(vecindad - 1, num_ciudades, seed, dist, solucion), 0, solucion, 0,
						solucion.length); // copiamos la solucion al meternos recursivamente
			} else {
				Random dado = new Random();
				int tirada = dado.nextInt(100);

				if (tirada < 50) { // si la tirada
					System.arraycopy(TSPLocalSearch(vecindad - 1, num_ciudades, seed, dist, solucion), 0, auxiliar, 0,
							auxiliar.length);

					if (evaluate(auxiliar, dist) < evaluate(solucion, dist)) {
						for (int j = 0; j < num_ciudades; j++) { // Para copiar un array a otro
							solucion[j] = auxiliar[j];
						}
					}
				} else { // Salimos del nivel de recursividad si falla la tirada
					fallo = true;
				}
			}
			vecindad--;
		}
		return solucion;
	}

	static void swap(int[] arr, int x, int y) {
		int temp = arr[x];
		arr[x] = arr[y];
		arr[y] = temp;
	}

	static void permute(int[] arr, double[][] dist, int[] solucion) {
		permute(arr, 0, arr.length - 1, dist, solucion);
	}

	static void permute(int[] arr, int i, int n, double[][] dist, int[] solucion) {
		int j;
		if (i == n) {
			if (evaluate(arr, dist) < evaluate(solucion, dist)) { // Si es mejor la guardamos
				for (int x = 0; x <= arr.length - 1; x++) {
					solucion[x] = arr[x];
				}
			}
		} else {
			for (j = i; j <= n; j++) {
				swap(arr, i, j);
				if (i != 0) {
					int[] arrayAuxiliar = rellenarAuxiliar(i, arr); // creamos un array auxiliar que va a tener desde 0
																	// hasta el elemento que hemos permutado, si al
																	// evaluar este array la distancia es mayor que la
																	// solucion podamos este camino
					if (evaluateAuxiliar(arrayAuxiliar, dist) < evaluate(solucion, dist)) { // Si el camino es
																							// prometedor seguimos si no
																							// no
						permute(arr, i + 1, n, dist, solucion); // nos adentramos un nivel mas, puede resultar que al
																// final sea mas largo pero en el paso anterior no lo
																// era
					}
					swap(arr, i, j); // deshacemos el swap
				} else { // Fallaba en el 0 porque el arrayauxiliar no podia coger el primer elemento
							// solo y decir si era menor, pero como siempre que sea principio de camino va a
							// ser menor dejamos que se meta un nivel de recursividad
					permute(arr, i + 1, n, dist, solucion);
					swap(arr, i, j); // backtrack
				}
			}
		}
	}

	/**
	 * Algoritmo Voráz
	 * 
	 * @param distances
	 * @param initialCity
	 * @return
	 * @throws Exception
	 */
	public static int[] GreedyAlgorithm(double[][] distances, int initialCity) throws Exception {// es la solucion
		// Genera un array solución para introducir los datos de nuestra solucion
		int[] route = new int[distances.length];
		// Numero de ciudades del problema
		int numCities = distances.length;

		// Indica por cual ciudad vamos a empezar
		int current = initialCity;

		int cont = 0;
		Stack<Integer> seleccionados = new Stack<Integer>();
		// int[] cities = new int[numCities];
		// Comprueba en la matriz de soluciones cuales son las mejores posiciones a las
		// que ir
		do {

			// if (current >= numCities) {
			// current =0;
			// for(int i=0;i<numCities;i++) {
			// if(route[current]==0) {
			// current=current+1;
			// }
			// }
			// }
			route[cont] = TSP_2.closestWay(distances, seleccionados, current, numCities);
			//
			// current++
			current=route[cont];
			cont++;

		} while (cont < numCities);

		for (int k = 0; k < route.length; k++) {
			System.out.println(route[k]);
		}
		return route;
	}

	public static int closestWay(double[][] distances, Stack<Integer> forbidden, int position, int ciudades)
			throws Exception {

		int res = 0;
		double pos = 0;

		double menor = Double.MAX_VALUE;
		for (int i = 0; i < ciudades; i++) {
			pos = distances[position][i];
			forbidden.push(i);
			if (!yaEscogida(forbidden, ciudades, position) && pos != 0 && pos < menor) {
				menor = pos;
				forbidden.push(i);
				res = i;
			} else {

				forbidden.pop();
			}

		}

		return res;

	}

	/**
	 * Verifica si la ciudad se encuentra dentro de las escogidas previamente
	 * 
	 * @param x
	 * @param numCities
	 * @param ciudad
	 * @return
	 * @throws Exception
	 */
	public static boolean yaEscogida(Stack<Integer> x, int numCities, int ciudad) throws Exception {
		
		boolean escogido = false;

		Iterator it = x.iterator();

		int actual = 0;
		while (it.hasNext()) {
			actual = (int) it.next();
			if (actual == ciudad) {
				escogido = true;
			}

		}

		return escogido;
	}

	public static int bestWay(double[][] b, Stack<Integer> x, int i) throws Exception {
		double min;
		int mejorTarea = -1;
		int position = 0;
		min = Double.MAX_VALUE;

		int noNulos = 0;
		// for(int l=0; l<x.length;l++) {
		// if(x[l]==0) {
		// noNulos++;
		// }
		// }

		// Nos aseguramos de que no vamos a tomar la direccion a si misma
		for (int tarea = 0; tarea < b.length; tarea++) {

			x.push(tarea);
			if (!TSP_2.yaEscogida(x, i, tarea) && b[i][tarea] < min && b[i][tarea] != 0) {
				min = b[i][tarea];
				mejorTarea = tarea;

				// x[position]=min;
				// position++;
			} else {
				x.pop();
			}
		}
		return mejorTarea;
	}

	/**
	 * 
	 * @param distances
	 * @param initialPosition
	 * @return
	 */
	public static int[] solve(double[][] distances, int initialPosition) {
		int[] route = new int[distances.length];
		int numCities = distances.length;

		int i = 1;
		int current = initialPosition;
		int[] citiesAvailable = new int[numCities];
		Arrays.fill(citiesAvailable, 0);
		route[0] = current;

		while (i < numCities) {
			citiesAvailable[current] = 1;
			int next = nextCity(citiesAvailable, current, distances, i);
			route[i] = next;
			current = next;
			i++;
		}
		return route;
	}

	private static int nextCity(int[] cavailable, int current, double[][] distances, int cRecorridas) {
		int nc = -1, j = 0;
		double path = Double.MAX_VALUE;
		double[][] dCiudad = new double[cavailable.length - cRecorridas][2], dOrdenado;

		for (int i = 0; i < cavailable.length; i++) {
			if (cavailable[i] == 0) {
				if (i < current) {

					dCiudad[j][0] = i;
					dCiudad[j][1] = distances[i][current - i];
					j++;

				} else if (i >= current) {
					dCiudad[j][0] = i;
					dCiudad[j][1] = distances[current][i - current];
					j++;
				}
			}
		}

		dOrdenado = divideYVenceras(dCiudad);

		nc = (int) dOrdenado[0][0];
		cavailable[nc] = 1;
		return nc;
	}

	// Funciones auxiliares
	public static double[][][] divide(double[][] entrada) {
		double[][][] s = new double[2][entrada.length / 2][2];
		int j = 0, k = 0;
		if (entrada.length % 2 == 1) {
			s[1] = new double[entrada.length / 2 + 1][2];
		}
		for (int i = 0; i < entrada.length; i++) {
			if (i < (entrada.length / 2)) {
				s[0][j] = entrada[i];
				j++;
			} else {
				s[1][k] = entrada[i];
				k++;
			}
		}
		return s;
	}

	public static double[][] combina(double[][][] entrada, int k) {
		double[][] s = new double[k][2];
		int j = 0, l = 0;
		for (int i = 0; i < entrada[0].length; i++) {
			while (entrada[1].length > j && entrada[0][i][1] > entrada[1][j][1]) {
				s[l] = entrada[1][j];
				j++;
				l++;
			}
			s[l] = entrada[0][i];
			l++;
		}
		if (l < k) {
			while (entrada[1].length > j) {
				s[l] = entrada[1][j];
				j++;
				l++;
			}
		}
		return s;
	}

	public static double pathCost(int[] path, double[][] dist) {
		double w = 0;

		for (int i = 0; i < path.length - 1; i++) {
			if (path[i + 1] < path[i]) {
				w += dist[path[i + 1]][path[i] - path[i + 1]];
			} else {
				w += dist[path[i]][path[i + 1] - path[i]];
			}
		}
		return w;
	}

	private static int nextCity(int[] cavailable, int current, double[][] distances) {
		int nc = -1;
		double path = Double.MAX_VALUE;

		for (int i = 0; i < cavailable.length; i++) {
			if (i < current) {
				if (cavailable[i] == 0 && distances[i][current - i] < path) {
					nc = i;
					path = distances[i][current - i];
				}
			} else if (i >= current) {
				if (cavailable[i] == 0 && distances[current][i - current] < path) {
					nc = i;
					path = distances[current][i - current];
				}
			}
		}
		cavailable[nc] = 1;
		return nc;
	}

	public static int[] rellenarAuxiliar(int j, int[] arr) {
		int[] aleatorio = creaMatriz(j);
		for (int i = 0; i < j; i++) {
			aleatorio[i] = arr[i];
		}
		return aleatorio;
	}

	public static double evaluate(int[] cities, double[][] dist) { // Esta funcion te va a dar el coste de la ruta
		double cost = 0.;

		for (int i = 0; i < cities.length - 1; i++) {
			int ii = cities[i];
			int jj = cities[i + 1];
			cost += dist[ii][jj];
		}

		cost += dist[cities.length - 1][0];

		return cost;
	}

	public static double evaluateAuxiliar(int[] cities, double[][] dist) { // Esta funcion te va a dar el coste de la
																			// ruta
		double cost = 0.;

		for (int i = 0; i < cities.length - 1; i++) {
			int ii = cities[i];
			int jj = cities[i + 1];
			cost += dist[ii][jj];
		}
		return cost;
	}

	public static int[] cambioCiudades(int[] recorrido, int n, int seed) { // Esta funcion va a generar un camino nuevo
																			// aleatorio generando dos numeros
																			// aleatorios y cambiando sus indices entre
																			// ellos

		Random rnd = new Random();
		Random rna = new Random();
		int indice1 = rnd.nextInt(n);
		int indice2 = rna.nextInt(n);
		int aux = 0;

		aux = recorrido[indice1];
		recorrido[indice1] = recorrido[indice2];
		recorrido[indice2] = aux;

		return recorrido;
	}

	public static int[] creaMatriz(int nc) { // Metodo para crear una matriz
		int[] A = new int[nc];

		for (int i = 0; i < A.length; i++) {
			A[i] = i;
		}
		return A;
	}

	public static double[] creaMatrizD(int nc) { // Metodo para crear una matriz
		double[] A = new double[nc];

		for (int i = 0; i < A.length; i++) {
			A[i] = (double) i;
		}
		return A;
	}

	private static boolean lineaContieneCoordenada(String lineaActual) {
		boolean equ = false;
		if (lineaActual.startsWith("1") || lineaActual.startsWith("2") || lineaActual.startsWith("3")
				|| lineaActual.startsWith("4") || lineaActual.startsWith("5") || lineaActual.startsWith("6")
				|| lineaActual.startsWith("7") || lineaActual.startsWith("8") || lineaActual.startsWith("9")) {
			equ = true;
		}
		return equ;
	}

	/**
	 * Método de lectura de ficheros TSP que devuelve la matriz de distancias del
	 * archivo TSP correspondiente.
	 * 
	 * @param FILENAME
	 *            nombre del archivo
	 * @return
	 */
	private static double[][] leerFichero(String FILENAME) {
		double[][] dist = null;
		double[][] ciudades = null;
		int nciudades = 0;
		BufferedReader br = null;
		FileReader fr = null;

		try {
			fr = new FileReader(FILENAME);
			br = new BufferedReader(fr);

			String lineaActual;

			br = new BufferedReader(new FileReader(FILENAME));
			lineaActual = br.readLine();
			while (lineaActual != null) {

				if (lineaActual.startsWith("DIMENSION")) {
					String[] lineaDividida2 = lineaActual.split(" ");
					nciudades = Integer.parseInt(lineaDividida2[lineaDividida2.length - 1]);
					ciudades = new double[nciudades][2];
					dist = new double[nciudades][nciudades];

				} else if (lineaActual.startsWith("NODE_COORD_SECTION")) {
					String[] lineaDividida2;
					int i = 0;
					lineaActual = br.readLine();
					while (lineaActual != null && i < nciudades) {
						int a = 0;
						lineaDividida2 = lineaActual.split(" ");

						while (!lineaContieneCoordenada(lineaDividida2[a])) {
							a++;
						}
						int ciudad = Integer.parseInt(lineaDividida2[a]);
						a++;

						while (!lineaContieneCoordenada(lineaDividida2[a])) {
							a++;
						}
						double x = Double.parseDouble(lineaDividida2[a]);
						a++;

						while (!lineaContieneCoordenada(lineaDividida2[a])) {
							a++;
						}
						double y = Double.parseDouble(lineaDividida2[a]);

						ciudades[ciudad - 1][0] = x;
						ciudades[ciudad - 1][1] = y;

						i++;
						lineaActual = br.readLine();
					}
				} else if (lineaActual != null) {
					System.out.println(lineaActual);
				}
				lineaActual = br.readLine();
			}

		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null) {
					br.close();
				}
				if (fr != null) {
					fr.close();
				}
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}

		for (int i = 0; i < nciudades; i++) {
			dist[i][i] = 0;
			for (int j = i + 1; j < nciudades; j++) {
				dist[i][j] = dist[j][i] = Math
						.sqrt((ciudades[i][0] - ciudades[j][0]) * (ciudades[i][0] - ciudades[j][0])
								+ (ciudades[i][1] - ciudades[j][1]) * (ciudades[i][1] - ciudades[j][1]));
			}
		}
		return dist;
	}
}
