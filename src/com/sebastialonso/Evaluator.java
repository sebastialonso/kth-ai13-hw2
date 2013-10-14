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

    /**
     * Performs the beta-pass algorithm
     * @return A Vector<Vector<Double>> with the rows being each beta_t
     */
   public Double[][] betaPass(){
        Double[][] betaMatrix = new Double[numberOfObservations][numberOfStates];

        for (int i=0; i < numberOfStates; i++){
            betaMatrix[numberOfObservations - 1][i] = 0.0;
        }

        for (int t = numberOfObservations-2; t >= 0; t--){
            for (int i=0; i < numberOfStates; i++){
                Double value = 0.0;
                for (int j=0; j < numberOfStates; j++){
                    value += Extended.esum( value,Extended.eproduct(
                                    Extended.eln(transitionMatrix[i][j]),
                                                 Extended.eproduct(
                                                    emissionMatrix[j][Integer.parseInt(observationsVector[t+1])],
                                                    Extended.eln(betaMatrix[t+1][j]))));
                }
                betaMatrix[t][i] = value;
            }
        }
        return betaMatrix;
   }
}


