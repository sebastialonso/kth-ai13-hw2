package com.sebastialonso;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Vector;

public class Main {

    public static void main(String[] args) throws Exception{
	// write your code here
        //FileInputStream stream = new FileInputStream("/home/seba/kth/ai13/hw2/kth-ai-hmm2-sample-data/hmm2_01.in");
        //FileInputStream stream = new FileInputStream("/home/seba/kth/ai13/hw2/kth-ai-hmm4-sample-data/hmm4_02.in");
        //System.setIn(stream);
        BufferedReader br = new BufferedReader(
                new InputStreamReader(System.in));



        Vector<String> lines = new Vector<String>();
        while (br.ready()){
            lines.add(br.readLine());
        }

        //Data structures from file
        Double[][] transitionMatrix = buildMatrix(lines.get(0));
        Double[][] emissionMatrix = buildMatrix(lines.get(1));
        Double[] initialVector = buildVector(lines.get(2));
        String[] observationVector = buildObservationVector(lines.get(3));

        //Evaluator evaluator = new Evaluator(transitionMatrix, emissionMatrix, initialVector, observationVector);
        //System.out.println(evaluator.evaluate());

        //Decoder decoder = new Decoder(transitionMatrix, emissionMatrix, initialVector, observationVector);
        //System.out.println(decoder.decode());
        //Evaluator evaluator = new Evaluator(transitionMatrix, emissionMatrix, initialVector, observationVector);
        //System.out.println(evaluator.evaluate());

        Learner learner = new Learner(transitionMatrix, emissionMatrix, initialVector, observationVector);
        //System.out.println(learner.learn(100));
        System.out.println(learner.learnReload(30));
    }

    /**
     * Builds a matrix out of a string with matrix information
     * @param matrixLine The String that contains number of rows, number of columns and elements
     * @return A double[][] as a matrix
     */
    private static Double[][] buildMatrix(String matrixLine){
        String[] matrixContent = matrixLine.split(" ");
        Double[][] matrix = new Double[Integer.parseInt(matrixContent[0])][Integer.parseInt(matrixContent[1])];
        int colIndex = 0;
        int rowIndex = 0;
        for (int i=2; i < matrixContent.length; i++){
            if (rowIndex == 0 || rowIndex%(Integer.parseInt(matrixContent[0])) != 0){
                matrix[rowIndex][colIndex] = Double.parseDouble(matrixContent[i]);
                colIndex++;
                if (colIndex != 0 && colIndex%(Integer.parseInt(matrixContent[1])) == 0){
                    rowIndex++;
                    colIndex = 0;
                }
            }
        }
        return matrix;
    }

    /**
     * Buils a vector with the observations
     * @param vectorLine The String that contains the number of observations and the observations
     * @return String[] that contains the observations as String elements
     */
    private static String[] buildObservationVector(String vectorLine){
        String[] vectorContent = vectorLine.split(" ");
        String[] observations = new String[Integer.parseInt(vectorContent[0])];

        for (int i=1; i< vectorContent.length; i++){
            observations[i-1] = vectorContent[i];
        }
        return observations;
    }

    /**
     * Buils a initial state vector for the HMM
     * @param vectorLine The String that contains the number of rows, columns and the probability distribution
     * @return A double[] with the initial probability distribution
     */
    private static Double[] buildVector(String vectorLine){
        String[] vectorContent = vectorLine.split(" ");
        Double[] vector = new Double[Integer.parseInt(vectorContent[1])];
        for (int i=2; i< vectorContent.length; i++ ){
            vector[i-2] = Double.parseDouble(vectorContent[i]);
        }
        return vector;
    }

    public static String matrixToString(Double[][] mat){
        String st= mat.length + " " + mat[0].length +" ";
        for (Double[] row : mat){
            for (Double element : row){
                st += element + " ";
            }
        }

        return st;
    }

    /**
     *
     * @param matrixes
     * @return
     */
    public static String printMatrixes(String[] matrixes){
        String st = "";
        for (String mat : matrixes){
            st += mat + "\n";
        }

        return st;
    }

    /**
     * DEBUG: Prints a matrix
     * @param matrix
     * @return
     */
    private static String printMatrix(Double[][] matrix){
        String st = "{\n";

        for (Double[] row : matrix){
            st += "{ ";
            for (Double element : row){
                st += element + " ";
            }
            st += "}\n";
        }
        st += " }";
        return st;

    }

    /**
     * DEBUG: Prints a vector
     * @param vector
     * @return
     */
    public static String printVector(Double[] vector){
        String st= "{ ";
        for (Double element : vector ){
            st += element + " ";
        }
        return st + "}";
    }
}
