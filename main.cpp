#include <iostream>
#include <math.h>
#include <complex>
#include <random>
#include <unordered_map>

/*
 *  Quantum Zeno paradox simulator
 *
 *  Date: 25.4.2021
 *
 *  Author: Blaz Rolih
 *          @blaz-r
 *         
 * */


/*
 *  Class for simulating quantum spin
 * */
class Spin {
private:
    std::complex<double> state[2];
    char symbol;    // |+> or |->
    double plusProbability;
    double minusProbability;
public:
    /*
     *  Constructor, creates wave function at given time
     *  At time 0 it's in state |+> and it's precessing due to field in Z direction
     *
     *  Equation: |Sig(t)> = 1/sqrt(2) * ( e^(-i*w*T/2) * up + e^(i*w*T/2) * down)
     *  where omega (w) is 2 * PI, up is [1, 0] and down is [0, 1]
     * */
    Spin(double time) {
        using namespace std::complex_literals;  // imaginary value (i)
        // makes state as given in equation above
        double exponent = 2 * M_PI * time / 2;
        std::complex<double> sqrt2 = 1 / sqrt(2);
        state[0] = sqrt2 * std::exp(-1i * exponent);
        state[1] = sqrt2 * std::exp(1i * exponent);
        symbol = '?';
    }
    /*
     *  Constructor that makes state from given x and y values
     * */
    Spin(double s1, double s2) { state[0] = s1; state[1] = s2; }

    /*
     *  Dot product operator for two spins, used for probability calculation
     * */
    std::complex<double> operator*(Spin spin2) {
        return this->state[0] * spin2.state[0] + this->state[1] * spin2.state[1];
    }

    /*
     *  Calculates probability that current spin is in state plus |+>
     * */
    std::complex<double> calculatePlusProbability() {
        // state plus |+> also 1/sqrt(2) * (up + down)
        Spin plus = Spin(1 / sqrt(2), 1 / sqrt(2));

        // Born rule, probability that current wave function is in state a: p(a) = |<a|Sig(t)>|^2
        std::complex<double> probability = plus * (*this);
        return abs(probability * probability);
    }

    /*
     *  Calculates probability that current spin is in state minus |->
     * */
    std::complex<double> calculateMinusProbability (){
        // state minus |-> also 1/sqrt(2) * (up - down)
        Spin minus = Spin(1 / sqrt(2), -1 / sqrt(2));

        // Born rule, probability that current wave function is in state a: p(a) = |<a|Sig(t)>|^2
        std::complex<double> probability = minus * (*this);
        return abs(probability * probability);
    }

    /*
     *  Measures current state of spin, which also collapses wave function
     * */
    char measure() {
        plusProbability = calculatePlusProbability().real();
        minusProbability = calculateMinusProbability().real();

        /*
         *  Generating random numbers according to discrete distribution,
         *  meaning that either |+> or |-> is selected depending on their probability
         * */
        std::random_device rnd;
        std::default_random_engine eng(rnd());
        std::discrete_distribution<> distr({plusProbability, minusProbability});
        double randomNum = distr(eng);

        // disc. distr. where index 0 is |+> and index 1 |->
        if(randomNum == 0) {
            symbol = '+';
            // collapse to |+>
            state[0] = 1 / sqrt(2);
            state[1] = 1 / sqrt(2);
        } else {
            symbol = '-';
            // collapse to |->
            state[0] = 1 / sqrt(2);
            state[1] = -1 / sqrt(2);
        }

        return symbol;
    }

    // symbol + / - getter
    char getSymbol() {
        return symbol;
    }

    // probability of + getter
    double getPlusProbability() {
        return plusProbability;
    }

    // probability of - getter
    double getMinusProbability() {
        return minusProbability;
    }
};

/*
 *  Class for demonstration of quantum Zeno paradox
 * */
class Zeno {
private:
    std::unordered_map<char, int> states;   // map containing number of each state

public:
    Zeno() {}

    /*
     *  Repeat the measurement of state selected number of time with given interval in between
     *  To print each measurement as + or -, print needs to be set on true
     * */
    void measure(bool print, double interval, int numberOfMeasurements) {
        // clear the map
        states['+'] = 0;
        states['-'] = 0;
        std::cout << "Measurements for time interval: " << interval << " periode." << std::endl;

        // start at 0 and increase time for given interval each time
        double time = 0;
        for (int i = 0; i < numberOfMeasurements; i++) {
            Spin waveFunction = Spin(time); // new wave function since it collapses
            char measurement = waveFunction.measure();  // + or - depending on time
            states[measurement]++;  // increase count of measured state
            if(print)
                std::cout << measurement << " ";
            time += interval;
        }
        if(print)
            std::cout << std::endl;

        double percentageOfPlus = states['+'] / (double) numberOfMeasurements * 100.0;
        double percentageOfMinus = 100 - percentageOfPlus;

        std::cout << "Number of plus measurements: " << states['+'] << ", " <<
                  percentageOfPlus << "%."<< std::endl;
        std::cout << "Number of minus measurements: " << states['-'] << ", " <<
                  percentageOfMinus << "%."<< std::endl;
    }

    /*
     *  Simulates quantum Zeno paradox by decreasing interval of measurement towards 0
     *  The closer it gets to 0, the less wave function changes,
     *  at the end we just keep measuring |+>
     * */
    void paradox(int numberOfMeasurements) {
        std::cout << "Quantum Zeno paradox: " << std::endl;
        // start at 1 period and decrease until 10^-10
        double interval = 1;
        while(interval > pow(10, -10)) {
            std::cout << "<-------------------------------------------->" << std::endl;
            measure(false, interval, numberOfMeasurements); // measure at given interval

            interval /= 2;
        }
    }
};

int main() {
    // measurement at 1 period interval
    Zeno* zeno1 = new Zeno();
    zeno1->measure(false, 1, 100000);
    std::cout << "<-------------------------------------------->" << std::endl;
    // measurement at 0.5 period interval
    Zeno* zeno05 = new Zeno();
    zeno05->measure(false, 0.5, 100000);
    std::cout << "<-------------------------------------------->" << std::endl;
    // measurement at 0.25 period interval
    Zeno* zeno025 = new Zeno();
    zeno025->measure(false, 0.25, 100000);
    std::cout << "<-------------------------------------------->" << std::endl;

    // simulation of quantum zeno paradox
    Zeno* zenoParadox = new Zeno();
    zenoParadox->paradox(100000);

    return 0;
}
