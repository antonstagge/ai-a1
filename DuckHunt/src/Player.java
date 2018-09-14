
class Player {


    final int TRAIN = 70;
    final double SHOOTING_THRESHOLD = 0.8;
    final double GUESSING_THRESHOLD = 0.8;
    final double STORK_THRESHOLD = -120;
    final int MAX_ITER = 500;
    final int NUMBER_OF_HIDDEN_STATES = 4;

    HMM[] species;

    HMM[] hmms;
    int timeStep = 0;
    int round = -1;
    int numberOfBirds = -1;

    public Player() {
        species = new HMM[Constants.COUNT_SPECIES];
    }

    /**
     * Shoot!
     *
     * This is the function where you start your work.
     *
     * You will receive a variable pState, which contains information about all
     * birds, both dead and alive. Each bird contains all past moves.
     *
     * The state also contains the scores for all players and the number of
     * time steps elapsed since the last time this function was called.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return the prediction of a bird we want to shoot at, or cDontShoot to pass
     */
    public Action shoot(GameState pState, Deadline pDue) {

        /*
         * Here you should write your clever algorithms to get the best action.
         * This skeleton never shoots.
         */

        if (round != pState.getRound()) {
            ++round;
            timeStep = 0;
            numberOfBirds = pState.getNumBirds();
            hmms = new HMM[numberOfBirds];
            for (int birdNumber = 0; birdNumber < numberOfBirds; ++birdNumber) {
                hmms[birdNumber] = new HMM(NUMBER_OF_HIDDEN_STATES, Constants.COUNT_MOVE);
            }
        }
        ++timeStep;


        if (timeStep < TRAIN) {
            return cDontShoot;
        }


        int[] bestMoves = new int[numberOfBirds];
        double[] bestMoveProbs = new double[numberOfBirds];

        for (int birdNumber = 0; birdNumber < numberOfBirds; ++birdNumber){
            Bird bird = pState.getBird(birdNumber);
            if (bird.isDead()) continue;

            if (timeStep == TRAIN) {
                hmms[birdNumber].trainHMM(MAX_ITER, getBirdFlightPath(pState.getBird(birdNumber)));
            } else {
                hmms[birdNumber].updateSequence(getBirdFlightPath(pState.getBird(birdNumber)));
            }

            //System.err.println("ONE HMM TRAINED:\n A for birdnumber: " + birdNumber + " is:");
            //System.err.println(hmms[birdNumber].A);

            double best_move_prob = -Double.MAX_VALUE;
            int best_move = -1;
            double[] pred = new double[Constants.COUNT_MOVE];
            for (int i = 0; i < Constants.COUNT_MOVE; ++i) {
                pred[i] = hmms[birdNumber].calcProbNextObs(i);
                if (pred[i] > best_move_prob) {
                    best_move_prob = pred[i];
                    best_move = i;
                }
            }

            bestMoveProbs[birdNumber] = best_move_prob;
            bestMoves[birdNumber] = best_move;

        }

        double bestProbOfAllBirds = -Double.MAX_VALUE;
        int bestMoveOfAllBirds = -1;
        int bestBird = -1;

        for (int birdNumber = 0; birdNumber < numberOfBirds; ++birdNumber) {
            Bird bird = pState.getBird(birdNumber);
            if (bird.isDead()) continue;
            if (bestMoveProbs[birdNumber] > bestProbOfAllBirds) {
                bestProbOfAllBirds = bestMoveProbs[birdNumber];
                bestMoveOfAllBirds = bestMoves[birdNumber];
                bestBird = birdNumber;
            }
        }

        if (bestProbOfAllBirds > SHOOTING_THRESHOLD) {
            double storkFitness = 0.0;
            if (species[Constants.SPECIES_BLACK_STORK] != null) {
                storkFitness = species[Constants.SPECIES_BLACK_STORK].speciesFitness(getBirdFlightPath(pState.getBird(bestBird)));

            }

            //System.err.println("STORK FITNESS: " + storkFitness);
            if (storkFitness < STORK_THRESHOLD) {
                Action shot = new Action(bestBird, bestMoveOfAllBirds);
                System.err.println("Shooting " + shot.toString() + " with prob: " + bestProbOfAllBirds);
                return shot;
            } else {
                //System.err.println("Stork danger!");
            }

        }


        // This line chooses not to shoot.
        //System.err.println("Ms left on return: " + pDue.remainingMs());
        return cDontShoot;

    }

    /**
     * Guess the species!
     * This function will be called at the end of each round, to give you
     * a chance to identify the species of the birds for extra points.
     *
     * Fill the vector with guesses for the all birds.
     * Use SPECIES_UNKNOWN to avoid guessing.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return a vector with guesses for all the birds
     */
    public int[] guess(GameState pState, Deadline pDue) {
        /*
         * Here you should write your clever algorithms to guess the species of
         * each bird. This skeleton makes no guesses, better safe than sorry!
         */

        
        if (round == 1) {
            int[] lGuess = new int[numberOfBirds];
            for (int i = 0; i < pState.getNumBirds(); ++i)
                lGuess[i] = Constants.SPECIES_PIGEON;
            return lGuess;
        }

        int[] lGuess = new int[numberOfBirds];
        for (int i = 0; i < numberOfBirds; ++i) {
            double maxFit = -Double.MAX_VALUE;
            int maxFitIdx = 0;
            for (int j = 0; j < Constants.COUNT_SPECIES; ++j) {
                if (species[j] == null) continue;
                int[] seq = getBirdFlightPath(pState.getBird(i));
                int[] new_seq = trimBirdFlightPath(seq);
                double fit = species[j].speciesFitness(new_seq);
                if (fit > maxFit) {
                    maxFit = fit;
                    maxFitIdx = j;
                }
            }
            lGuess[i] = maxFitIdx;

        }

        return lGuess;
    }

    /**
     * If you hit the bird you were trying to shoot, you will be notified
     * through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pBird the bird you hit
     * @param pDue time before which we must have returned
     */
    public void hit(GameState pState, int pBird, Deadline pDue) {
        System.err.println("HIT BIRD!!!");
    }

    /**
     * If you made any guesses, you will find out the true species of those
     * birds through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pSpecies the vector with species
     * @param pDue time before which we must have returned
     */
    public void reveal(GameState pState, int[] pSpecies, Deadline pDue) {
        for (int i = 0; i < pState.getNumBirds(); ++i) {
            if (species[pSpecies[i]] == null)
                if (pSpecies[i] == Constants.SPECIES_BLACK_STORK) {
                    System.err.println("                STORK DETECTED");
                }
                species[pSpecies[i]] = new HMM(hmms[i]);

            int[] seq = getBirdFlightPath(pState.getBird(i));
            int[] new_seq = trimBirdFlightPath(seq);
            species[pSpecies[i]].trainHMM(MAX_ITER, new_seq);
        }
    }

    public static final Action cDontShoot = new Action(-1, -1);

    public int[] getBirdFlightPath(Bird b) {
        int[] seq = new int[b.getSeqLength()];
        for (int i = 0; i < b.getSeqLength(); ++i) {
            seq[i] = b.getObservation(i);
        }
        return seq;
    }

    public int[] trimBirdFlightPath(int[] seq) {
        int l = 0;
        for (int k = 0; k < seq.length; ++k) {
            if (seq[k] == -1)
                break;
            ++l;
        }
        int[] new_seq = new int[l];
        for (int k = 0; k < l; ++k) {
            new_seq[k] = seq[k];
        }
        return new_seq;
    }
}
