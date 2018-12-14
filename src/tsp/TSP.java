package tsp;

import java.io.*;
import java.util.*;


public class TSP {
	//Version ultima 
	public static int[] generaPrimera(int numCiudades) {
        int[] sol = new int[numCiudades];
        for (int i = 0; i < numCiudades; i++) {
            sol[i] = i;
        }
        return sol;
    }
	
	public static double fuerzaBruta(int[] ciudades, int i, int n, double mejor, double[][] dist) {
        int j;
        double aux;
        if (i == n) {
            aux = getCoste(ciudades, dist);
            if (aux < mejor) {
                System.out.println("");
                for (int k = 0; k < ciudades.length; k++) {
                    System.out.print(ciudades[k] + ",");
                    if (k % 10 == 0 && k != 0) {
                        System.out.println("");
                    }
                }
                mejor = aux;
            }
        } else {
            for (j = i; j < n; j++) {
                swap(ciudades, i, j);
                mejor = fuerzaBruta(ciudades, i + 1, n, mejor, dist);
            }
        }
        return mejor;
    }
	
	 public static int[] swap(int[] array, int x, int y) {
	        int aux = array[x];
	        array[x] = array[y];
	        array[y] = aux;
	        return array;
	    }
    
    public static double getCoste(int[] ciudad, double[][] distancia) {
    	double coste = 0;
    	
    	for (int i = 0; i < ciudad.length-1; i++) {
			int ii=ciudad[i];
			int jj=ciudad[i+1];
			coste+=distancia[ii][jj];
		}
    	
    	coste+=distancia[ciudad.length-1][0];
    	
    	return coste;
    }
    
    public static int[] generaSolucion(int n) {
        int[] s1 = new int[n];
        for (int i = 0; i < n; i++) {
            s1[i] = i;
        }
        Random rand = new Random();

        for (int i = 0; i < n; i++) {
            int a = rand.nextInt(n);
            int b = rand.nextInt(n);
            int aux = s1[a];
            s1[a] = s1[b];
            s1[b] = aux;
        }
        return s1;
    }

    public static int[] busquedaAleatoriaConIteracionesFijas(int nCiudades, double[][] distancia, int busquedas) {

        int[] s1 = new int[nCiudades];
        int[] s2 = new int[nCiudades];
        for (int i = 0; i < nCiudades; i++) {
            s1[i] = i;
        }

        for (int i = 0; i < busquedas; i++) {
            s2 = generaSolucion(nCiudades);
            double a = getCoste(s1, distancia);
            double b = getCoste(s2, distancia);
            if (b < a) {
                for (int j = 0; j < nCiudades; j++) {
                    s1[j] = s2[j];
                }
            }
        }
        return s1;
    }

    public static int[] busquedaAleatoriaIteracionesReiniciables(int nCiudades, double[][] distancia, int busquedas) {
        int[] s1 = new int[nCiudades];
        int[] s2 = new int[nCiudades];
        for (int i = 0; i < nCiudades; i++) {
            s1[i] = i;
        }
        int b = 0;
        while (b < busquedas) {
            s2 = generaSolucion(nCiudades);
            if (getCoste(s2, distancia) < getCoste(s1, distancia)) {
                for (int i = 0; i < nCiudades; i++) {
                    s1[i] = s2[i];
                }
                b = 0;
            } else {
                b++;
            }
        }
        return s1;
    }
    
    private static double[][] leerFichero(String nombreFichero) {
        double[][] distancia = null;
        double[][] ciudades = null;
        int nCiudades = 0;
        BufferedReader br = null;
        FileReader fr = null;

        try {
            fr = new FileReader(nombreFichero);
            br = new BufferedReader(fr);

            String lineaActual;

            br = new BufferedReader(new FileReader(nombreFichero));
            lineaActual = br.readLine();
            while (lineaActual != null) {

                if (lineaActual.startsWith("DIMENSION")) {
                    String[] lineaDividida2 = lineaActual.split(" ");
                    nCiudades = Integer.parseInt(lineaDividida2[lineaDividida2.length - 1]);
                    System.out.println(lineaActual);
                    ciudades = new double[nCiudades][2];
                    distancia = new double[nCiudades][nCiudades];

                } else if (lineaActual.startsWith("NODE_COORD_SECTION")) {
                    String[] lineaDividida2;
                    int i = 0;
                    lineaActual = br.readLine();
                    while (lineaActual != null && i < nCiudades) {
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
        
        for (int i = 0; i < nCiudades; i++) {
            distancia[i][i] = 0;
            for (int j = i + 1; j < nCiudades; j++) {
                distancia[i][j] = distancia[j][i] = Math.sqrt(
                        (ciudades[i][0] - ciudades[j][0]) * (ciudades[i][0] - ciudades[j][0])
                        + (ciudades[i][1] - ciudades[j][1]) * (ciudades[i][1] - ciudades[j][1]));
            }
        }
        return distancia;
    }
    
    
    private static boolean lineaContieneCoordenada(String lineaActual) {
        boolean equ = false;
        if (lineaActual.startsWith("1") || lineaActual.startsWith("2")
                || lineaActual.startsWith("3") || lineaActual.startsWith("4")
                || lineaActual.startsWith("5") || lineaActual.startsWith("6")
                || lineaActual.startsWith("7") || lineaActual.startsWith("8")
                || lineaActual.startsWith("9")) {
            equ = true;
        }
        return equ;
    }

 
    public static void main(String[] args) {
    	
    	System.out.println("\nEXPERIMENTO DEL TSP\n\n");
        int iteraciones = 100;
        String[] instancias = {"berlin52.tsp", "kroA100.tsp", "kroA150.tsp", "kroA200.tsp",
                "a280.tsp", "vm1084.tsp", "vm1748.tsp", "usa13509.tsp"};

        double[][] duraciones = new double[instancias.length][3];

        for (int i = 0; i < instancias.length; i++) {
            System.out.println("\n\n--------------------------------------------");
            System.out.println("Instancia " + (i + 1));
            double[][] distancia = leerFichero("data\\" + instancias[i]);
            int nCiudades = distancia.length;
            double s = 0;
            int[] s1 = new int[nCiudades];
            int repeticiones = 50;
            

            System.out.println("\nEjecutando Fuerza bruta");
            double start = System.currentTimeMillis();
            	s = fuerzaBruta(s1, 0, 9, getCoste(s1, distancia), distancia);
            double stop = System.currentTimeMillis();
            duraciones[i][0] = (stop - start);

            System.out.println("Ejecutando Búsqueda aleatoria iteraciones fijas");
            start = System.currentTimeMillis();
            for (int j = 0; j < repeticiones; j++) {
                s1 = busquedaAleatoriaConIteracionesFijas(nCiudades, distancia, iteraciones);
            }
            stop = System.currentTimeMillis();
            duraciones[i][1] = (stop - start) / repeticiones;

            System.out.println("Ejecutando Búsqueda aleatoria iteraciones reiniciables");
            start = System.currentTimeMillis();
            for (int j = 0; j < repeticiones; j++) {
                s1 = busquedaAleatoriaIteracionesReiniciables(nCiudades, distancia, iteraciones);
            }
            stop = System.currentTimeMillis();
            duraciones[i][2] = (stop - start) / repeticiones;
            System.out.println("--------------------------------------------");
            
        }

        System.out.println("\n\n\n--------------------------------------------");
        System.out.println("Resultados finales");
        for (int i1 = 0; i1 < instancias.length; i1++) {
            System.out.println("Instancia " + (i1 + 1));
            System.out.println("\tFB: " + duraciones[i1][0]+" ms");     //Fuerza bruta
            System.out.println("\tBAIF: " + duraciones[i1][1]+" ms");  //Búsqueda aleatoria con iteraciones fijas
            System.out.println("\tBAIR: " + duraciones[i1][2]+" ms");  //Búsqueda aleatoria con iteraciones reiniciables
        }
    

        
    }
} 
    
    
    