package com.sebastialonso;

import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: Sebastián González Mardones
 * Date: 10/8/13
 * Time: 10:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class Evaluator {

    private Vector<Vector<Double>> transitionMatrix;
    private Vector<Vector<Double>> emissionMatrix;
    private Vector<Double> initialVector;
    private Vector<String>  observationsVector;


    public Evaluator(Vector<Vector<Double>> transition, Vector<Vector<Double>> emission, Vector<Double> initial, Vector<String> observations){
        this.transitionMatrix = transition;
        this.emissionMatrix = emission;
        this.initialVector = initial;
        this.observationsVector = observations;
    }

    public Double evaluate(){
        return sumElements(alphaPass().lastElement());
    }

    public Vector<Vector<Double>> alphaPass(){
        int numberOfStates = this.transitionMatrix.size();
        int numberOfTimes = this.observationsVector.size();

        Vector<Vector<Double>> alphaMatrix = new Vector<Vector<Double>>(numberOfTimes);
        Vector<Double> alphaZero = new Vector<Double>();
        for (int j=0; j < numberOfStates; j++){

            alphaZero.add(this.initialVector.get(j) * this.emissionMatrix.get(j).get(Integer.parseInt(this.observationsVector.elementAt(0))));
        }
        alphaMatrix.add(alphaZero);
        Vector<Vector<Double>> transitionTranspose = transpose(this.transitionMatrix);
        for (int t=1; t < numberOfTimes; t++){
            String currentObservation = this.observationsVector.elementAt(t);
            alphaMatrix.add(alpha_t(alphaMatrix.get(t-1), numberOfStates, transitionTranspose, this.emissionMatrix,currentObservation));
        }
        return alphaMatrix;

    }
    /**
     * Performs once cycle for a time t, updating the value of the alpha array
     * @param previousAlpha Vector<Double> of the previous time
     * @param numberOfStates Number of slots in prevoiusAlpha (number of states of the HMM)
     * @param transpose The transpose of the transition matrix, used for dot product
     * @param emission The emission matrix, B.
     * @param currentObservation The String corresponding to the observation O_t
     * @return
     */
    public static Vector<Double> alpha_t(Vector<Double> previousAlpha, int numberOfStates,
                                          Vector<Vector<Double>> transpose, Vector<Vector<Double>> emission,
                                          String currentObservation ){
        Vector<Double> newAlpha = new Vector<Double>();
        for (Integer i=0; i< numberOfStates; i++){
            newAlpha.add(dot(previousAlpha, transpose.get(i)) * emission.get(i).get(Integer.parseInt(currentObservation)));
        }
        return newAlpha;
    }

    /**
     * Performs the beta-pass algorithm
     * @return A Vector<Vector<Double>> with the rows being each beta_t
     */
    public Vector<Vector<Double>> betaPass(){
        int numberOfObservations =  this.observationsVector.size();
        int numberOfStates = this.transitionMatrix.size();
        //Initialize beta
        Vector<Vector<Double>> betaMatrix = new Vector<Vector<Double>>();
        Vector<Double> betaZero = new Vector<Double>(numberOfStates);
        for (int index=0; index < numberOfStates; index++){
            betaZero.add(1.0);
        }
        betaMatrix.add(betaZero);
        for (int t = numberOfObservations-2; t >= 0; t--){
            String futureObservation = this.observationsVector.elementAt(t+1);
            System.out.println("t = " + t + " | " + "O_" + (t + 1) + ": " + futureObservation);
            betaMatrix.insertElementAt(beta_t(betaMatrix.get(betaMatrix.size() - 1), numberOfStates, this.transitionMatrix, this.emissionMatrix, futureObservation), 0);
        }
        return betaMatrix;
    }

    private static Vector<Double> beta_t(Vector<Double> posteriorBeta, int numberOfStates,
                                         Vector<Vector<Double>> transition, Vector<Vector<Double>> emission, String futureObservation){
        Vector<Double> currentBeta = new Vector<Double>();
        for (int i=0; i < numberOfStates; i++){
            for (int j=0; j < numberOfStates; j++){
                currentBeta.add(transition.get(i).get(j) * emission.get(j).get(Integer.parseInt(futureObservation)) * posteriorBeta.get(j));
            }
        }
        return currentBeta;
    }

    /**
     * Transposes a matrix, i.e, changes its rows per its columns
     * @param matrix A Vector<Vector<Double>> matrix to be transposed
     * @return A Vector<Vector<Double>> matrix which has rows as the columns of A, and columns as the rows of A
     */
    public static Vector<Vector<Double>> transpose(Vector<Vector<Double>> matrix){
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

    /**
     * Sums the element of a vector
     * @param vector Vector<Double> on which the internal sum is desired
     * @return Double sum of the elements of the vector
     */
    private static Double sumElements(Vector<Double> vector){
        Double sum = 0.0;
        for (Double element : vector){
            sum += element;
        }

        return sum;
    }

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
}

