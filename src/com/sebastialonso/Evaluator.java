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

    private double[][] transitionMatrix;
    private double[][] emissionMatrix;
    private double[] initialVector;
    private double[] scalingFactor;
    private String[]  observationsVector;
    private int numberOfStates;
    private int numberOfObservations;


    public Evaluator(double[][] transition, double[][] emission, double[] initial, String[] observations){
        this.transitionMatrix = transition;
        this.emissionMatrix = emission;
        this.initialVector = initial;
        this.observationsVector = observations;
        this.numberOfStates = transitionMatrix.length;
        this.numberOfObservations = observationsVector.length;
    }

    public Double evaluate(){
        return alphaPass();
    }

    /**
     * Used to solve Problem 2: Evaluation. No scaling needed
     * @return
     */
    public Double alphaPass(){
        double[][] alpha = new double[numberOfObservations][numberOfStates];

        scalingFactor[0] = 0.0;

        for (int i=0; i < numberOfStates; i++){
            alpha[0][i] = initialVector[i] * emissionMatrix[i][Integer.parseInt(observationsVector[0])];
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
                    alpha[t][i] += alpha[t-1][j] * transitionMatrix[j][i];
                }
                alpha[t][i] *= emissionMatrix[i][Integer.parseInt(observationsVector[t])];
                scalingFactor[t] += alpha[t][i];
            }
            //Scale a_t(i)
            scalingFactor[t] = 1/scalingFactor[t];
            for (int i=0; i < numberOfStates; i++){
                alpha[t][i] *= scalingFactor[t];
            }
        }
         double prob = 0;
        for (int i=0; i < numberOfObservations; i++){
            prob += Math.log(scalingFactor[i]);
        }
        return  -prob;
    }

    /**
     * Used to solve Problem 4: Learn. Scaling needed
     * @param scalingFactor Vector<Double> where the scaling factors are stored.
     * @return A Matrix with the alpha Vector for each T
     */
    /*public Vector<Vector<Double>> alphaPass(Vector<Double> scalingFactor){
        Vector<Vector<Double>> alphaMatrix = new Vector<Vector<Double>>();

        Vector<Double> alphaZero = new Vector<Double>();
        Double scale = 0.0;
        for (int i=0; i < numberOfStates; i++){

            alphaZero.add(initialVector.get(i) * emissionMatrix.get(i).get(Integer.parseInt(observationsVector.get(0))));
            scale += initialVector.get(i) * emissionMatrix.get(i).get(Integer.parseInt(observationsVector.get(0)));
        }

        //Scaling the vector
        scalingFactor.add( 1/ scale);
        for (int i=0; i < numberOfStates;i++){
            alphaZero.set(i, alphaZero.get(i) * scalingFactor.get(0));
        }

        alphaMatrix.add(alphaZero);

        for (int t=1; t < numberOfObservations; t++){
            int currentObservation = Integer.parseInt(observationsVector.get(t));
            Vector<Double> newAlpha= new Vector<Double>();
            scale = 0.0;
            for (int i=0; i< numberOfStates; i++){
                Double value = 0.0;
                for (int j=0; j < numberOfStates; j++){
                    value +=  alphaMatrix.get(t-1).get(j) * transitionMatrix.get(j).get(i);
                }
                value *= emissionMatrix.get(i).get(currentObservation);
                scale += value;
                newAlpha.add(value);
            }
            //Scaling the whole vector
            scalingFactor.add(1/ scale);
            for (int i =0; i < numberOfStates; i++){
                newAlpha.set(i, newAlpha.get(i) * scalingFactor.get(t));
            }

            alphaMatrix.add(newAlpha);
        }
        return alphaMatrix;
    } */
    /**
     * Performs the beta-pass algorithm
     * @return A Vector<Vector<Double>> with the rows being each beta_t
     */
   /* public Vector<Vector<Double>> betaPass(Vector<Double> scalingFactor){
        Vector<Vector<Double>> betaMatrix = new Vector<Vector<Double>>(numberOfObservations);

        Vector<Double> betaZero = new Vector<Double>(numberOfStates);
        for (int i=0; i < numberOfStates; i++){
            betaZero.add(scalingFactor.lastElement());
        }


        betaMatrix.add(betaZero);

        for (int t = numberOfObservations-2; t >= 0; t--){
            String futureObservation = observationsVector.elementAt(t+1);
            Vector<Double> currentBeta = new Vector<Double>();
            for (int i=0; i < numberOfStates; i++){
                Double value = 0.0;
                for (int j=0; j < numberOfStates; j++){
                    value += transitionMatrix.get(i).get(j) * emissionMatrix.get(j).get(Integer.parseInt(futureObservation)) * betaMatrix.get(betaMatrix.size() - 1).get(j);
                }
                currentBeta.add(value);

                //Scale b_t
                currentBeta.set(i, currentBeta.get(i) * scalingFactor.get(t));
            }
            betaMatrix.add(betaMatrix.size() - 1 ,currentBeta);
        }
        return betaMatrix;
    }*/

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

