package com.sebastialonso;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.Vector;

public class Main {

    public static void main(String[] args) throws Exception{
	// write your code here
        //FileInputStream stream = new FileInputStream("/home/seba/kth/ai13/hw2/kth-ai-hmm2-sample-data/hmm2_01.in");
        FileInputStream stream = new FileInputStream("/home/seba/kth/ai13/hw2/kth-ai-hmm4-sample-data/hmm4_01.in");
        System.setIn(stream);
        BufferedReader br = new BufferedReader(
                new InputStreamReader(System.in));



        Vector<String> lines = new Vector<String>();
        while (br.ready()){
            lines.add(br.readLine());
        }

        //Data structures from file
        Vector<Vector<Double>> transitionMatrix = buildMatrix(lines.get(0));
        Vector<Vector<Double>> emissionMatrix = buildMatrix(lines.get(1));
        Vector<Double> initialVector = buildVector(lines.get(2));

        if (lines.size() == 4){
            Vector<String> observationVector = buildObservationVector(lines.get(3));

            //Decoder decoder = new Decoder(transitionMatrix, emissionMatrix, initialVector, observationVector);
            //System.out.println(decoder.decode());
            //Evaluator evaluator = new Evaluator(transitionMatrix, emissionMatrix, initialVector, observationVector);
            //System.out.println(evaluator.evaluate());

            Learner learner = new Learner(transitionMatrix, emissionMatrix, initialVector, observationVector);
            System.out.println(learner.learn(60));
            System.out.println(learner.learnReload(40));

        }


    }

    /**
     * Builds a matrix out of a string with matrix information
     * @param matrixLine The String that contains number of rows, number of columns and elements
     * @return A Vector<Vector<Double>> as a matrix
     */
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

    /**
     * Buils a vector with the observations
     * @param vectorLine The String that contains the number of observations and the observations
     * @return A Vector<String> that contains the observations as String elements
     */
    private static Vector<String> buildObservationVector(String vectorLine){
        String[] vectorContent = vectorLine.split(" ");
        Vector<String> observations = new Vector<String>(Integer.parseInt(vectorContent[0]));

        for (int i=1; i< vectorContent.length; i++){
            observations.add(vectorContent[i]);
        }
        return observations;
    }

    /**
     * Buils a initial state vector for the HMM
     * @param vectorLine The String that contains the number of rows, columns and the probability distribution
     * @return A Vector<Double> with the probability distribution
     */
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
     *//*
    private static Vector<Double> multiplyVectorMatrix(Vector<Double> vector, Vector<Vector<Double>> matrix){
        //Initialized with matrix col number
        Vector<Double> result = new Vector<Double>(matrix.get(0).size());
        Vector<Vector<Double>> transpose = transpose(matrix);

        for (int columnOfA = 0; columnOfA < matrix.get(0).size(); columnOfA++){
            result.add(dot(vector, transpose.get(columnOfA)));
        }
        return result;
    }

    *//**
     * Transposes a matrix, i.e, changes its rows per its columns
     * @param matrix A Vector<Vector<Double>> matrix to be transposed
     * @return A Vector<Vector<Double>> matrix which has rows as the columns of A, and columns as the rows of A
     *//*
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
    }*/

    /**
     * Calculates the dot product between two vectors
     * @param first A Vector<Double> to be multiplied
     * @param second The second Vector<Double> to be multiplied
     * @return A Double with the value of the dot product
     */
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

    public static String printVector(Vector<Double> vector){
        String st= "{ ";
        for (Double element : vector ){
            st += element + " ";
        }
        return st + "}";
    }

    private static String finalPrintVector(Vector<Double> vector){
        String st="1 " + vector.size();
        for (Double element : vector){
            st += " " + new DecimalFormat("#.##").format(element);
        }
        return st;

    }
}
