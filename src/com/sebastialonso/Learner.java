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


    public Learner(Vector<Vector<Double>> transition, Vector<Vector<Double>> emission, Vector<Double> initial, Vector<String> observations){
        this.transitionMatrix = transition;
        this.emissionMatrix = emission;
        this.initialState = initial;
        this.observationsVector = observations;
        this.numberOfStates = transition.size();
        this.numberOfObservations = observations.size();
    }

    public Object[] estimateModel(){
        Evaluator evaluator = new Evaluator(transitionMatrix, emissionMatrix, initialState, observationsVector);

        Vector<Vector<Double>> alpha = evaluator.alphaPass();
        Vector<Vector<Double>> beta = evaluator.betaPass();
        Vector<Vector<Vector<Double>>> xiMatrix = createXi(alpha,beta);
        Vector<Vector<Double>> gammaMatrix = createGamma(xiMatrix);


        return null;
    }

    private Vector<Vector<Double>> createGamma(Vector<Vector<Vector<Double>>> xiMatrix){
        Vector<Vector<Double>> gammaMatrix = new Vector<Vector<Double>>();
        for (int t=0; t < numberOfObservations; t++){
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
        for (int t=0; t < numberOfObservations; t++){
            eachXi = xiMatrix(alpha, beta, t);
            hyperMatrix.add(eachXi);
        }
        return hyperMatrix;
    }

    private Double denominator(Vector<Vector<Double>> alpha, Vector<Vector<Double>> beta, int currentT){
        Double sum = 0.0;
        for (int i=0; i < numberOfStates; i++){
            for (int j=0; j < numberOfStates; i++){
                sum += alpha.get(currentT).get(i) * transitionMatrix.get(i).get(j) * emissionMatrix.get(j).get(Integer.parseInt(this.observationsVector.elementAt(currentT+1))) * beta.get(currentT + 1).get(j);
            }
        }
        return sum;
    }
}
