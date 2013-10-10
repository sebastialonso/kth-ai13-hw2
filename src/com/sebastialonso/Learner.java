package com.sebastialonso;

import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: Sebastián González Mardones
 * Date: 10/8/13
 * Time: 10:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class Learner {
    private Vector<Vector<Double>> transitionMatrix;
    private Vector<Vector<Double>> emissionMatrix;
    private Vector<Double> initialState;
    private Vector<String> observationsVector;
    private int numberOfStates;
    private int numberOfObservations;
    private int numberOfSymbols;


    public Learner(Vector<Vector<Double>> transition, Vector<Vector<Double>> emission, Vector<Double> initial, Vector<String> observations){
        this.transitionMatrix = transition;
        this.emissionMatrix = emission;
        this.initialState = initial;
        this.observationsVector = observations;
        this.numberOfStates = transition.size();
        this.numberOfObservations = observations.size();
        this.numberOfSymbols = emission.get(0).size();
    }


    public String learn(int numberOfIterations){
        //Get lambda_0
        Double DELTA = 1e-8;

        Evaluator previousModel = new Evaluator(this.transitionMatrix, this.emissionMatrix, this.initialState, this.observationsVector);
        Vector<Vector<Double>> alpha = previousModel.alphaPass();
        Vector<Vector<Double>> beta = previousModel.betaPass();

        for (int iteration=0; iteration < numberOfIterations; iteration++){

            Vector<Vector<Double>> transitionEstimation = estimateTransition(alpha, beta);
            Vector<Vector<Double>> emissionEstimation = estimateEmission(alpha, beta);
            Vector<Double> initialVectorEstimation = estimatePi(alpha, beta);
            Evaluator nextModel = new Evaluator(transitionEstimation, emissionEstimation, initialVectorEstimation, this.observationsVector);

            //End Condition
            Double previousEstimation = previousModel.evaluate();
            Double currentEstimation = nextModel.evaluate();
            if (Math.abs( previousEstimation - currentEstimation) < DELTA){
                System.out.println();
                Vector<Vector<Vector<Double>>> answer = new Vector<Vector<Vector<Double>>>();
                answer.add(transitionEstimation);
                answer.add(emissionEstimation);
                return printMatrixes(answer);
            }
            else {
                alpha = nextModel.alphaPass();
                beta = nextModel.betaPass();
                previousModel = nextModel;
                System.out.println("Previous Estimation: " + previousEstimation);
                System.out.println("Current Estimation: " + currentEstimation);
            }
        }
        //Run out of computations


        return  "Iterations completed. No convergence found";
    }

    /**
     * Estimates the Emission Matrix B
     * @param alpha
     * @param beta
     * @return
     */
    private Vector<Vector<Double>> estimateEmission(Vector<Vector<Double>> alpha, Vector<Vector<Double>> beta){

        Vector<Vector<Vector<Double>>> xiMatrix = createXi(alpha,beta);
        Vector<Vector<Double>> gammaMatrix = createGamma(xiMatrix);

        Vector<Vector<Double>> emission = new Vector<Vector<Double>>();
        for (int j=0; j < numberOfStates; j++){
            Vector<Double> emissionRow = new Vector<Double>();
            for (int k=0; k < numberOfSymbols; k++){
                Double numerator = 0.0;
                Double denominator = 0.0;
                for (int t=0; t < numberOfObservations -2; t++){
                    denominator += gammaMatrix.get(t).get(j);
                    if (Integer.parseInt(observationsVector.get(t)) == k){
                        numerator += gammaMatrix.get(t).get(j);
                    }
                }
                emissionRow.add(numerator/denominator);
            }
            emission.add(emissionRow);
        }
        return emission;
    }

    /**
     * Estimates the Transition Matrix A
     * @param alpha
     * @param beta
     * @return
     */
    private Vector<Vector<Double>> estimateTransition(Vector<Vector<Double>> alpha, Vector<Vector<Double>> beta){

        Vector<Vector<Vector<Double>>> xiMatrix = createXi(alpha,beta);
        Vector<Vector<Double>> gammaMatrix = createGamma(xiMatrix);

        Vector<Vector<Double>> transition = new Vector<Vector<Double>>();
        for (int i=0; i < numberOfStates; i++){
            Vector<Double> estimationVector = new Vector<Double>();
            for (int j=0; j < numberOfStates; j++){
                Double numerator = 0.0;
                Double denominator = 0.0;
                for (int t=0; t < numberOfObservations - 2; t++){
                    numerator += xiMatrix.get(t).get(i).get(j);
                    denominator += gammaMatrix.get(t).get(i);

                }
                estimationVector.add(numerator/denominator);
            }
            transition.add(estimationVector);
        }

        return transition;

    }

    /**
     * Estimates the Initial Vector Pi
     * @param alpha
     * @param beta
     * @return
     */
    private Vector<Double> estimatePi(Vector<Vector<Double>> alpha, Vector<Vector<Double>> beta){
        Vector<Vector<Vector<Double>>> xiMatrix = createXi(alpha,beta);
        Vector<Vector<Double>> gammaMatrix = createGamma(xiMatrix);

        Vector<Double> pi = new Vector<Double>();
        for (int i=0; i < gammaMatrix.get(0).size(); i++){
            pi.add(gammaMatrix.get(0).get(i));
        }
        return pi;
    }

    /**
     *
     * @param xiMatrix
     * @return
     */
    private Vector<Vector<Double>> createGamma(Vector<Vector<Vector<Double>>> xiMatrix){
        Vector<Vector<Double>> gammaMatrix = new Vector<Vector<Double>>();
        for (int t=0; t < numberOfObservations - 2; t++){
            Vector<Double> gamma = new Vector<Double>();
            for (int i=0; i < numberOfStates; i++){
                Double valueAtI = 0.0;
                for (int j=0; j < numberOfStates; j++){
                    valueAtI += xiMatrix.get(t).get(i).get(j);
                }
                gamma.add(valueAtI);
            }
            gammaMatrix.add(gamma);

        }
        return gammaMatrix;
    }

    /**
     * Creates a xi matrix given a time t (as defined in hmmmut.pdf page 26)
     * @param alpha Matrix from alpha-pass
     * @param beta Matrix from beta-pass
     * @return Vector<Vector<Double>> Xi matrix
     */
    private Vector<Vector<Double>> xiMatrix(Vector<Vector<Double>> alpha, Vector<Vector<Double>> beta, int currentT){
        Vector<Vector<Double>> xiMatrix = new Vector<Vector<Double>>();
        Vector<Double> xi = new Vector<Double>();

        for (int i=0; i < numberOfStates; i++){
            for (int j=0; j < numberOfStates; j++){
                Double numerator = alpha.get(currentT).get(i) * transitionMatrix.get(i).get(j) *
                        emissionMatrix.get(j).get(Integer.parseInt(observationsVector.elementAt(currentT+1))) *
                        beta.get(currentT+1).get(j);
                xi.add(numerator / denominator(alpha,beta,currentT));
            }
            xiMatrix.add(xi);
        }

        return xiMatrix;
    }

    /**
     * Creates a Vector of Xi Matrices, one matrix per time of observation
     * @param alpha
     * @param beta
     * @return
     */
    private Vector<Vector<Vector<Double>>> createXi(Vector<Vector<Double>> alpha, Vector<Vector<Double>> beta){
        Vector<Vector<Vector<Double>>> hyperMatrix = new Vector<Vector<Vector<Double>>>();
        Vector<Vector<Double>> eachXi;
        for (int t=0; t < numberOfObservations - 2; t++){
            eachXi = xiMatrix(alpha, beta, t);
            hyperMatrix.add(eachXi);
        }
        return hyperMatrix;
    }

    private Double denominator(Vector<Vector<Double>> alpha, Vector<Vector<Double>> beta, int currentT){
        Double sum = 0.0;
        for (int i=0; i < numberOfStates; i++){
            for (int j=0; j < numberOfStates; j++){
                sum += alpha.get(currentT).get(i) * transitionMatrix.get(i).get(j) * emissionMatrix.get(j).get(Integer.parseInt(this.observationsVector.elementAt(currentT+1))) * beta.get(currentT + 1).get(j);
            }
        }
        return sum;
    }

    private String printMatrixes(Vector<Vector<Vector<Double>>> matrixes){
        String st = "";
        for (Vector<Vector<Double>> matrix : matrixes){
            st += matrix.size() + " " + matrix.get(0).size() + " ";
            for (Vector<Double> row : matrix){
                for (Double element : row){
                    st += element + " ";
                }
            }
            st += "\n";
        }

        return st;

    }
}

