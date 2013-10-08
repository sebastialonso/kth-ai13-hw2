package com.sebastialonso;

import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: Sebastián González Mardones
 * Date: 10/8/13
 * Time: 10:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class Decoder {

    private Vector<Vector<Double>> transitioMatrix;
    private Vector<Vector<Double>> emissionMatrix;
    private Vector<Double> initialVector;
    private Vector<String>  observationsVector;

    public Decoder( Vector<Vector<Double>> transition, Vector<Vector<Double>> emission, Vector<Double> initial, Vector<String> observations){
        this.transitioMatrix = transition;
        this.emissionMatrix = emission;
        this.initialVector = initial;
        this.observationsVector = observations;
    }

    public String decode(){
        int numberOfStates = this.transitioMatrix.size();
        int numberOfObservations = this.observationsVector.size();

        Vector<Vector<Double>> deltaMatrix = new Vector<Vector<Double>>();
        Vector<Double> delta = new Vector<Double>();
        Vector<Vector<Integer>> phiMatrix = new Vector<Vector<Integer>>();
        Vector<Integer> phi = new Vector<Integer>();

        for (int i=0; i < numberOfStates; i++){
            delta.add(this.initialVector.get(i) * this.emissionMatrix.get(i).get(Integer.parseInt(this.observationsVector.elementAt(0))));
            phi.add(0);
        }
        deltaMatrix.add(delta);
        phiMatrix.add(phi);

        //Recursion part
        for (int t=1; t < numberOfObservations; t++){
            Vector<Double> newDelta = new Vector<Double>();
            Vector<Integer> newPhi = new Vector<Integer>();
            String currentObservation = this.observationsVector.elementAt(t);

            for (int j=0; j< numberOfStates; j++){
                Double maxValue = -1.0;
                Integer index = -1;
                for (int i=0; i < numberOfStates; i++){
                    Double value = deltaMatrix.get(t-1).get(i) * this.transitioMatrix.get(i).get(j);
                    if (value > maxValue){
                        maxValue = value;
                        index = i;
                    }
                }
                newDelta.add(maxValue * this.emissionMatrix.get(j).get(Integer.parseInt(currentObservation)));
                newPhi.add(index);
            }
            deltaMatrix.add(newDelta);
            phiMatrix.add(newPhi);


        }

        //Getting the indexes
        Integer index = -1;
        Double maxValue = -1.0;
        for (int i=0; i < numberOfStates; i++){
            if (deltaMatrix.lastElement().get(i) > maxValue){
                maxValue = deltaMatrix.lastElement().get(i);
                index = i;
            }
        }
        int[] states = new int[numberOfObservations];
        states[numberOfObservations - 1] = index;
        for (int t = numberOfObservations - 2; t >= 0; t--){
            states[t] = phiMatrix.get(t + 1).get(states[t + 1]);
        }

        //Printing the indexes
        String st = "";
        for (int i= 0; i < numberOfObservations; i++){
            st += states[i] + " ";
        }

        return st;
    }
}
