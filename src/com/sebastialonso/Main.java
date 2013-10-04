package com.sebastialonso;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Vector;

public class Main {

    public static void main(String[] args) throws Exception{
	// write your code here
        BufferedReader br = new BufferedReader(
                new InputStreamReader(System.in));
        String[] line = new String[3];

        for (int i=0; i<3 ; i++){
            line[i] = br.readLine();
        }

        Vector<Vector<Double>> transitionMatrix = buildMatrix(line[0]);
        Vector<Vector<Double>> emissionMatrix = buildMatrix(line[1]);
        Vector<Double> initialVector = buildVector(line[2]);
        System.out.println(printMatrix(transitionMatrix));
        System.out.println(printVector(initialVector));
        System.out.println(printMatrix(transpose(transitionMatrix)));
        System.out.println("pi * A: \n" + multiplyVectorMatrix(initialVector, transitionMatrix));

    }

    private static Vector<Vector<Double>> buildMatrix(String matrixLine){
        String[] matrixContent = matrixLine.split(" ");
        Vector<Vector<Double>> matrix = new Vector<Vector<Double>>();
        int colIndex = 0;
        int rowIndex = 0;
        Vector<Double> row = new Vector<Double>();
        for (int i=2; i < matrixContent.length; i++){
            if (rowIndex == 0 || rowIndex%(Integer.parseInt(matrixContent[0])) != 0){
                //matrix.get(rowIndex).add(colIndex,Double.parseDouble(matrixContent[i]));
                row.add(Double.parseDouble(matrixContent[i]));
                colIndex++;
                if (colIndex != 0 && colIndex%(Integer.parseInt(matrixContent[1])) == 0){
                    rowIndex++;
                    colIndex = 0;
                    matrix.add(row);
                    row = new Vector<Double>();
                }
            }
        }
        return matrix;
    }

    private static Vector<Double> buildVector(String vectorLine){
        String[] vectorContent = vectorLine.split(" ");
        Vector<Double> vector = new Vector<Double>();
        for (int i=2; i< vectorContent.length; i++ ){
            vector.add(Double.parseDouble(vectorContent[i]));
        }
        return vector;
    }

    /**
     * Performs matrix multiplication, and returns the result
     * @param vector Vector that multiplies matrix
     * @param matrix Matrix being multiplied
     * @return result, vector with as many rows as vector and as many columns as matrix
     */
    private static Vector<Double> multiplyVectorMatrix(Vector<Double> vector, Vector<Vector<Double>> matrix){
        //Initialized with matrix col number
        Vector<Double> result = new Vector<Double>(matrix.get(0).size());
        Vector<Vector<Double>> transpose = transpose(matrix);

        for (int columnOfA = 0; columnOfA < matrix.get(0).size(); columnOfA++){
            result.add(dot(vector, transpose.get(columnOfA)));
        }
        return result;
    }

    /**
     * Transposes a matrix, i.e, changes its rows per its columns
     * @param matrix
     * @return
     */
    private static Vector<Vector<Double>> transpose(Vector<Vector<Double>> matrix){
        int rowNumber = matrix.size();
        int colNumber = matrix.get(0).size();
        Vector<Vector<Double>> transposeMatrix = new Vector<Vector<Double>>(rowNumber);
        //We scanned the whole row
        for (int i=0; i< colNumber;i++){
            Vector<Double> newRows = new Vector<Double>(colNumber);
            for (Vector<Double> row : matrix){
                newRows.add(row.get(i));
            }
            transposeMatrix.add(newRows);
        }
        return transposeMatrix;
    }

    private static Double dot(Vector<Double> first, Vector<Double> second){
        Double result = 0.0;
        for (int i=0; i< first.size(); i++){
            result += first.get(i) * second.get(i);
        }
        return result;
    }

    private static String printMatrix(Vector<Vector<Double>> matrix){
        String st = "{\n";

        for (Vector<Double> row : matrix){
            st += "{ ";
            for (Double element : row){
                st += element + " ";
            }
            st += "}\n";
        }
        st += " }";
        return st;

    }

    private static String printVector(Vector<Double> vector){
        String st= "{ ";
        for (Double element : vector ){
            st += element + " ";
        }
        return st + "}";
    }
}
