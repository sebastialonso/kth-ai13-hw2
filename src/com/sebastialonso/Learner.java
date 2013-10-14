package com.sebastialonso;


public class Learner {
    private Double[][] transitionMatrix;
    private Double[][] emissionMatrix;
    private Double[] initialState;
    private String[] observationsVector;
    private int numberOfStates;
    private int numberOfObservations;
    private int numberOfSymbols;
    private  double DELTA = 1e-13;


    public Learner(Double[][] transition, Double[][] emission, Double[] initial, String[] observations){
        this.transitionMatrix = transition;
        this.emissionMatrix = emission;
        this.initialState = initial;
        this.observationsVector = observations;
        this.numberOfStates = transition.length;
        this.numberOfObservations = observations.length;
        this.numberOfSymbols = emission[0].length;
    }

    /**public String learnReload(int iterations){
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
            double[][] createGamma = new double[numberOfObservations][numberOfStates];
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

            ///compute createGamma och diggama
            for (int t=0; t < numberOfObservations - 1; t++){
                double denominator = 0.0;
                for (int i=0; i < numberOfStates; i++){
                    for (int j=0; j< numberOfStates; j++){
                        denominator += alpha[t][i] * transition.get(i).get(j) * emission.get(j).get(Integer.parseInt(observationsVector.get(t+1))) * beta[t+1][j];
                    }
                }
                for (int i=0; i < numberOfStates; i++){
                    createGamma[t][i] = 0.0;
                    for (int j=0; j < numberOfStates; j++){
                        diGamma[t][i][j] = (alpha[t][i] * transition.get(i).get(j)*emission.get(j).get(Integer.parseInt(observationsVector.get(t+1))) * beta[t+1][j])/ denominator;
                        createGamma[t][i] += diGamma[t][i][j];
                    }
                }
            }

            ///Re-estimate model
            //Re.estimate pi
            for (int i=0; i < numberOfStates; i++){
                estimatedInitial.add(createGamma[0][i]);
            }

            //Re-estimate A
            for (int i=0; i < numberOfStates; i++){
                Vector<Double> estimatedRow = new Vector<Double>();
                for (int j=0; j < numberOfStates; j++){
                    double numerator = 0.0;
                    double denominator = 0.0;

                    for (int t=0; t < numberOfObservations-1; t++){
                        numerator += diGamma[t][i][j];
                        denominator += createGamma[t][i];
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
                            numerator += createGamma[t][i];
                        }
                        denominator += createGamma[t][i];
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
    }  */


    public String learn(int iterations){

        Double currentProb = Double.NEGATIVE_INFINITY;

        Evaluator currentModel = new Evaluator(transitionMatrix, emissionMatrix, initialState, observationsVector);
        Double[][] alpha = currentModel.alphaPass();
        Double[][] beta = currentModel.betaPass();
        String[] answer = new String[2];

        for (int iter = 0; iter < iterations; iter++){

            Double[][][] xi = createXi(alpha, beta);
            Double[][] gamma = createGamma(alpha, beta);


            Double[] estimatedInitial = estimatePi(gamma);
            Double[][] estimatedTransition = estimateTransition(xi, gamma);
            Double[][] estimatedEmission = estimateEmission(gamma);

            Evaluator newModel = new Evaluator(estimatedTransition, estimatedEmission, estimatedInitial, observationsVector);
            Double newProb = newModel.evaluate();

            if (Math.abs(currentProb - newProb) < DELTA){
                answer[0] = Main.matrixToString(estimatedTransition);
                answer[1] = Main.matrixToString(estimatedEmission);
                break;
            }
            else {
                currentProb = newProb;
                alpha = newModel.alphaPass();
                beta = newModel.betaPass();
            }
        }

        return Main.printMatrixes(answer);
    }

    /**
     * Compute Pi, the initial probability vector
     * @param gamma
     * @return
     */
    private Double[] estimatePi(Double[][] gamma){
        Double[] pi = new Double[numberOfStates];
        for (int i = 0; i < numberOfStates; i++){
            pi[i]  = Extended.eexp(gamma[0][i]);
        }

        return pi;
    }

    /**
     * Computes the estimation for the transition matrix
     * @param xi
     * @param gamma
     * @return
     */
    private Double[][] estimateTransition(Double[][][] xi, Double[][] gamma){

        Double[][] transition = new Double[numberOfStates][numberOfStates];
        Double numerator = Double.NaN;
        Double denominator = Double.NaN;

        for (int t = 0; t < numberOfObservations - 1; t++){
            for (int i = 0; i < numberOfStates; i++){
                for (int j = 0; j < numberOfStates; j++){
                    numerator = Extended.esum(numerator, xi[t][i][j]);
                    denominator = Extended.esum(denominator, gamma[t][i]);
                }
            }

        }

        for (int i = 0; i < numberOfStates; i++){
            for (int j = 0; j < numberOfStates; j++){
                transition[i][j] = Extended.eexp(
                        Extended.eproduct(numerator, -denominator)
                );
            }
        }

        return transition;
    }

    /**
     * Computes the estimation for the emission matrix
     * @param gamma
     * @return
     */
    private Double[][] estimateEmission(Double[][] gamma){
        Double[][] emission = new Double[numberOfStates][numberOfSymbols];

        Double numerator = Double.NaN;
        Double denominator = Double.NaN;

        for (int j = 0; j < numberOfStates; j++){
            for (int k = 0; k < numberOfSymbols; k++){
                for (int t = 0; t < numberOfObservations; t++){
                    if (Integer.parseInt(observationsVector[t]) == k){
                        numerator = Extended.esum(numerator, gamma[t][j]);
                    }
                    denominator = Extended.esum(denominator, gamma[t][j]);
                }
                emission[j][k] = Extended.eexp(Extended.eproduct(numerator, - denominator));
            }
        }

        return emission;
    }

    /**
     * Computes the log space Gamma Matrix
     * @param alpha
     * @param beta
     * @return
     */
    private Double[][] createGamma(Double[][] alpha, Double[][] beta){

        Double[][] gamma = new Double[numberOfObservations][numberOfStates];
        for (int t = 0; t < numberOfObservations; t++){
            Double normalizer = Double.NaN;
            for (int i = 0; i < numberOfStates; i++){
                gamma[t][i] = Extended.eproduct(alpha[t][i], beta[t][i]);
                normalizer = Extended.esum(normalizer, gamma[t][i]);
            }
            for (int i = 1; i < numberOfStates; i++){
                gamma[t][i] = Extended.eproduct(gamma[t][i], - normalizer);
            }
        }

        return gamma;
    }

    /**
     * Computes the log space Xi Matrix
     * @param alpha
     * @param beta
     * @return
     */
    private Double[][][] createXi(Double[][] alpha, Double[][] beta){
        Double[][][] xiMatrix = new Double[numberOfObservations][numberOfStates][numberOfStates];

        for (int t = 0; t < numberOfObservations - 1; t++){
            Double normalize = Double.NaN;
            for (int i = 0; i < numberOfStates; i++){
                for (int j = 0; j < numberOfStates; j++){
                    xiMatrix[t][i][j] = Extended.eproduct(
                            alpha[t][i],
                            Extended.eproduct(
                                    Extended.eln(transitionMatrix[i][j]),
                                    Extended.eproduct(
                                            Extended.eln(emissionMatrix[j][Integer.parseInt(observationsVector[t+1])]),
                                            beta[t+1][j]
                                    )
                            )
                    );
                    normalize = Extended.esum(normalize, xiMatrix[t][i][j]);
                }
            }

            for (int i = 0; i < numberOfStates; i++){
                for (int j = 0; j < numberOfStates; j++){
                    xiMatrix[t][i][j] = Extended.eproduct(xiMatrix[t][i][j], -normalize);
                }
            }
        }

        return xiMatrix;
    }
}

