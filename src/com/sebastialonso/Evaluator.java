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

    private Double[][] transitionMatrix;
    private Double[][] emissionMatrix;
    private Double[] initialVector;
    private String[]  observationsVector;
    private int numberOfStates;
    private int numberOfObservations;


    public Evaluator(Double[][] transition, Double[][] emission, Double[] initial, String[] observations){
        this.transitionMatrix = transition;
        this.emissionMatrix = emission;
        this.initialVector = initial;
        this.observationsVector = observations;
        this.numberOfStates = transitionMatrix.length;
        this.numberOfObservations = observationsVector.length;
    }

    public Double evaluate(){
        return sumElements(alphaPass()[numberOfObservations -1]);
    }

    /**
     * Used to solve Problem 2: Evaluation.
     * @return
     */
    public Double[][] alphaPass(){
        Double[][] alpha = new Double[numberOfObservations][numberOfStates];

        for (int i=0; i < numberOfStates; i++){
            alpha[0][i] = Extended.eproduct(
                    Extended.eln(initialVector[i]),
                    Extended.eln(emissionMatrix[i][Integer.parseInt(observationsVector[0])]));

        }

        //compute a_t(i)
        for (int t = 1; t< numberOfObservations; t++){
            for (int i=0; i < numberOfStates; i++){
                Double val = Double.NaN;
                for (int j=0; j < numberOfStates; j++){
                    val = Extended.esum(val, Extended.eproduct(alpha[t - 1][j], Extended.eln(transitionMatrix[j][i])));
                }
                alpha[t][i] = Extended.eproduct(val, Extended.eln(emissionMatrix[i][Integer.parseInt(observationsVector[t])]));

            }
        }

        return alpha;
    }

    /**
     * Sums the element of a vector
     * @param vector Vector<Double> on which the internal esum is desired
     * @return Double esum of the elements of the vector
     */
    public Double sumElements(Double[] vector){
        Double response = 0.0;
        for (int i=0; i< numberOfStates; i++){
            response = Extended.esum(response, vector[i]);
        }
        return Extended.eexp(response)-1;
    }
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


