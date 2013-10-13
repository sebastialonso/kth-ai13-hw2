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
    private  double DELTA = 1e-13;


    public Learner(Vector<Vector<Double>> transition, Vector<Vector<Double>> emission, Vector<Double> initial, Vector<String> observations){
        this.transitionMatrix = transition;
        this.emissionMatrix = emission;
        this.initialState = initial;
        this.observationsVector = observations;
        this.numberOfStates = transition.size();
        this.numberOfObservations = observations.size();
        this.numberOfSymbols = emission.get(0).size();
    }
    public String learnReload(int iterations){
        Double oldLogProb = Double.NEGATIVE_INFINITY;
        Vector<Vector<Double>> transition = transitionMatrix;
        Vector<Vector<Double>> emission = emissionMatrix;
        Vector<Double> initial = initialState;


        for (int iteration=0; iteration < iterations; iteration++){

            Vector<Double> estimatedInitial = new Vector<Double>();
            Vector<Vector<Double>> estimatedTransition = new Vector<Vector<Double>>();
            Vector<Vector<Double>> estimatedEmission = new Vector<Vector<Double>>();

            ///The alpha pass
            double[] scalingFactor = new double[numberOfObservations];
            double[][] alpha = new double[numberOfObservations][numberOfStates];
            double[][] beta = new double[numberOfObservations][numberOfStates];
            double[][] gamma = new double[numberOfObservations][numberOfStates];
            double[][][] diGamma = new double[numberOfObservations][numberOfStates][numberOfStates];


            //Compute alpha[0]

            scalingFactor[0] = 0.0;

            for (int i=0; i < numberOfStates; i++){
                alpha[0][i] = initial.get(i) * emission.get(i).get(Integer.parseInt(observationsVector.get(0)));
                scalingFactor[0] += alpha[0][i];
            }

            //Scale alpha[0]
            scalingFactor[0] = 1/scalingFactor[0];
            for (int i=0; i < numberOfStates; i++){
                alpha[0][i] *= scalingFactor[0];
            }

            //compute a_t(i)
            for (int t = 1; t< numberOfObservations; t++){
                scalingFactor[t] = 0.0;
                for (int i=0; i < numberOfStates; i++){
                    alpha[t][i] = 0.0;
                    for (int j=0; j < numberOfStates; j++){
                        alpha[t][i] += alpha[t-1][j] * transition.get(j).get(i);
                    }
                    alpha[t][i] *= emission.get(i).get(Integer.parseInt(observationsVector.get(t)));
                    scalingFactor[t] += alpha[t][i];
                }
                //Scale a_t(i)
                scalingFactor[t] = 1/scalingFactor[t];
                for (int i=0; i < numberOfStates; i++){
                    alpha[t][i] *= scalingFactor[t];
                }
            }

            ///The beta pass
            //Scale beta
            for (int i=0; i < numberOfStates; i++){
                beta[numberOfObservations - 1][i] = scalingFactor[numberOfObservations-1];
            }

            //beta pass
            for (int t= numberOfObservations - 2; t >= 0; t--){
                for (int i=0; i < numberOfStates; i++){
                    beta[t][i] = 0.0;
                    for (int j=0; j < numberOfStates; j++){
                        beta[t][i] = beta[t][i] + transition.get(i).get(j) * emission.get(j).get(Integer.parseInt(observationsVector.get(t+1))) * beta[t+1][j];
                    }
                    //scale beta_t
                    beta[t][i] *= scalingFactor[t];
                }
            }

            ///compute gamma och diggama
            for (int t=0; t < numberOfObservations - 1; t++){
                double denominator = 0.0;
                for (int i=0; i < numberOfStates; i++){
                    for (int j=0; j< numberOfStates; j++){
                        denominator += alpha[t][i] * transition.get(i).get(j) * emission.get(j).get(Integer.parseInt(observationsVector.get(t+1))) * beta[t+1][j];
                    }
                }
                for (int i=0; i < numberOfStates; i++){
                    gamma[t][i] = 0.0;
                    for (int j=0; j < numberOfStates; j++){
                        diGamma[t][i][j] = (alpha[t][i] * transition.get(i).get(j)*emission.get(j).get(Integer.parseInt(observationsVector.get(t+1))) * beta[t+1][j])/ denominator;
                        gamma[t][i] += diGamma[t][i][j];
                    }
                }
            }

            ///Re-estimate model
            //Re.estimate pi
            for (int i=0; i < numberOfStates; i++){
                estimatedInitial.add(gamma[0][i]);
            }

            //Re-estimate A
            for (int i=0; i < numberOfStates; i++){
                Vector<Double> estimatedRow = new Vector<Double>();
                for (int j=0; j < numberOfStates; j++){
                    double numerator = 0.0;
                    double denominator = 0.0;

                    for (int t=0; t < numberOfObservations-1; t++){
                        numerator += diGamma[t][i][j];
                        denominator += gamma[t][i];
                    }
                    estimatedRow.add(numerator/denominator);
                }
                estimatedTransition.add(estimatedRow);
            }

            //Re-estimate B
            for (int i=0; i < numberOfStates; i++){
                Vector<Double> estimatedEmissionRow = new Vector<Double>();
                for (int j=0; j < numberOfSymbols; j++){
                    double numerator = 0.0;
                    double denominator = 0.0;
                    for (int t=0; t < numberOfObservations -1; t++){
                        if (Integer.parseInt(observationsVector.get(t)) ==  j){
                            numerator += gamma[t][i];
                        }
                        denominator += gamma[t][i];
                    }
                    estimatedEmissionRow.add(numerator/denominator);
                }
                estimatedEmission.add(estimatedEmissionRow);
            }

            ///Compute log[P(O|lambda)]
            double logProb =0;
            for (int t=0; t < numberOfObservations; t++){
                logProb += Math.log(scalingFactor[t]);
            }
            logProb = -1 * logProb;

            //Move the values
            transition = estimatedTransition;
            emission = estimatedEmission;
            initial = estimatedInitial;

            if (iteration < iterations && Math.abs(logProb - oldLogProb) < DELTA){
                break;

            }
            else {
                oldLogProb = logProb;
            }



        }
        Vector<Vector<Vector<Double>>> response = new Vector<Vector<Vector<Double>>>();
        response.add(transition);
        response.add(emission);

        return printMatrixes(response);
    }

    public String learn(int iterations){
        Double previousProb = Double.NEGATIVE_INFINITY;

        Vector<Double> estimatedInitial = new Vector<Double>();
        Vector<Vector<Double>> estimatedTransition = new Vector<Vector<Double>>();
        Vector<Vector<Double>> estimatedEmission = new Vector<Vector<Double>>();

        Evaluator currentModel = new Evaluator(transitionMatrix, emissionMatrix, initialState, observationsVector);
        Vector<Double> scalingFactor = new Vector<Double>();
        Vector<Vector<Double>> alpha = currentModel.alphaPass(scalingFactor);
        Vector<Vector<Double>> beta = currentModel.betaPass(scalingFactor);

        Vector<Vector<Vector<Double>>> answer = new Vector<Vector<Vector<Double>>>();
        for (int iter=0; iter < iterations; iter++){


            Vector<Vector<Vector<Double>>> xi = createXi(alpha, beta);
            Vector<Vector<Double>> gamma = createGamma(xi);

            estimatedInitial = estimatePi(gamma);
            estimatedTransition = estimateTransition(xi, gamma);
            estimatedEmission = estimateEmission(gamma);

            Evaluator newModel = new Evaluator(estimatedTransition, estimatedEmission, estimatedInitial, observationsVector);

            Double newModelProb = newModel.evaluate();
            if (newModelProb > previousProb){
                answer.add(estimatedTransition);
                answer.add(estimatedEmission);
                break;
            }
            else{
                previousProb = newModelProb;
                scalingFactor = new Vector<Double>();
                alpha = newModel.alphaPass(scalingFactor);
                beta = newModel.betaPass(scalingFactor);
            }
        }

        return printMatrixes(answer);
    }

    /*public String learn1(int numberOfIterations){
        Double oldLogProb = Double.NEGATIVE_INFINITY;

        Evaluator previousModel = new Evaluator(this.transitionMatrix, this.emissionMatrix, this.initialState, this.observationsVector);

        Vector<Double> scalingFactor = new Vector<Double>();
        Vector<Vector<Double>> alpha = previousModel.alphaPass(scalingFactor);
        Vector<Vector<Double>> beta = previousModel.betaPass(scalingFactor);

        for (int iteration=0; iteration < numberOfIterations; iteration++){

            Vector<Vector<Double>> transitionEstimation = estimateTransition(alpha, beta);
            Vector<Vector<Double>> emissionEstimation = estimateEmission(alpha, beta);
            Vector<Double> initialVectorEstimation = estimatePi(alpha, beta);
            Evaluator nextModel = new Evaluator(transitionEstimation, emissionEstimation, initialVectorEstimation, this.observationsVector);

            //End Condition
            Double previousEstimation = previousModel.evaluate();
            Double currentEstimation = nextModel.evaluate();
            if (Math.abs( previousEstimation - currentEstimation) < DELTA){
                //System.out.println();
                Vector<Vector<Vector<Double>>> answer = new Vector<Vector<Vector<Double>>>();
                answer.add(transitionEstimation);
                answer.add(emissionEstimation);
                return printMatrixes(answer);
            }
            else {
                alpha = nextModel.alphaPass(scalingFactor);
                beta = nextModel.betaPass(scalingFactor);
                previousModel = nextModel;
                //System.out.println("Previous Estimation: " + previousEstimation);
                //System.out.println("Current Estimation: " + currentEstimation);
            }
        }
        return "asdads";
    }*/

    /**
     *
     * @param gammaMatrix
     * @return
     */
    private Vector<Vector<Double>> estimateEmission(Vector<Vector<Double>> gammaMatrix){

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
     *
     * @param xiMatrix
     * @param gammaMatrix
     * @return
     */
    private Vector<Vector<Double>> estimateTransition(Vector<Vector<Vector<Double>>> xiMatrix, Vector<Vector<Double>> gammaMatrix){

        Vector<Vector<Double>> transition = new Vector<Vector<Double>>();
        for (int i=0; i < numberOfStates; i++){
            Vector<Double> transitionRow = new Vector<Double>();
            Double numerator, denominator;
            for (int j=0; j < numberOfStates; j++){
                numerator = 0.0;
                denominator = 0.0;
                for (int t=0; t < numberOfObservations - 2; t++){
                    numerator += xiMatrix.get(t).get(i).get(j);
                    denominator += gammaMatrix.get(t).get(i);
                }
                transitionRow.add(numerator / denominator);
            }
            transition.add(transitionRow);
        }

        return transition;

    }

    /**
     *
     * @param gammaMatrix
     * @return
     */
    private Vector<Double> estimatePi(Vector<Vector<Double>> gammaMatrix){

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
            xi = new Vector<Double>();
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
        for (int t=0; t <= numberOfObservations - 2; t++){
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

    private Vector<Vector<Double>> normalizeAndSetScalingVector( Vector<Vector<Double>> matrix, Vector<Double> scalingFactor){
        Vector<Vector<Double>> normalized = new Vector<Vector<Double>>();
        Double sum = 0.0;
        for (Vector<Double> alphaRow : matrix){
            Vector<Double> newAlphaRow = new Vector<Double>();

            for (Double elem : alphaRow){
                sum += elem;
            }
            scalingFactor.add(1 / sum);
            for (Double elem : alphaRow){
                newAlphaRow.add( elem * 1/sum);
            }
            normalized.add(newAlphaRow);
        }

        return normalized;
    }
}

