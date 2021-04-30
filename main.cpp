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
 *  Requires c++ version 14
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
     *  Makes evolution for given time, ie. rotates state according to original equation
     *
     *  Equation: |Sig(t)> = 1/sqrt(2) * ( e^(-i*w*T/2) * up + e^(i*w*T/2) * down)
     *  where omega (w) is 2 * PI, up is [1, 0] and down is [0, 1]
     * */
    void evolution(double time) {
        using namespace std::complex_literals;  // imaginary value (i)
        // makes evolution for given time
        double exponent = 2 * M_PI * time / 2;
        state[0] *= std::exp(-1i * exponent);
        state[1] *= std::exp(1i * exponent);
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
        return conj(this->state[0]) * spin2.state[0] + conj(this->state[1]) * spin2.state[1];
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
    std::complex<double> calculateMinusProbability() {
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
    void measure(bool print, bool collapse, double interval, int numberOfMeasurements) {
        // clear the map
        states['+'] = 0;
        states['-'] = 0;
        std::cout << "Measurements for time interval: " << interval << " period." << std::endl;

        // start at 0 and increase time for given interval each time
        double time = 0;
        Spin waveFunction = Spin(time);
        for (int i = 0; i < numberOfMeasurements; i++) {
            if(collapse) {
                /* Since we measure the function every time, it collapses in either |+> or |->
                 * and we then make evolution from that collapsed state forwards.
                 * This is how it works in real life
                 * */
                waveFunction.evolution(time);
            } else {
                /* This is example for a universe, unlike ours, where there's no collapse of wave function
                 * We don't have any collapses, so we always just start in state [+> and then rotate it
                 * for given time, all according to equation defined in class Spin
                 * */
                waveFunction = Spin(time);
            }
            char measurement = waveFunction.measure();  // + or - depending on time
            states[measurement]++;  // increase count of measured state
            if(print)
                std::cout << time << "T: " << measurement << std::endl;
            time += interval;
        }

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
        // start at 1 period and decrease until 10^-12
        double interval = 1;
        while(interval > pow(10, -12)) {
            std::cout << "<-------------------------------------------->" << std::endl;
            measure(false, true, interval, numberOfMeasurements); // measure at given interval

            interval /= 2;
        }
    }
};

void testSameObject() {
    std::cout << "Test using evolution of SAME object:" << std::endl;
    std::cout << "First with time 1 period, then 0,5 then 0,25 and finally 0." << std::endl;

    Spin sameObj = Spin(1);
    std::cout << "[T = 1]: " << "p|+> = " << sameObj.calculatePlusProbability().real()
              << ", p|-> = " << sameObj.calculateMinusProbability().real() << std::endl;
    sameObj.evolution(0.5);
    std::cout << "[T = 0.5]: " << "p|+> = " << sameObj.calculatePlusProbability().real()
              << ", p|-> = " << sameObj.calculateMinusProbability().real() << std::endl;
    sameObj.evolution(0.25);
    std::cout << "[T = 0.25]: " << "p|+> = " << sameObj.calculatePlusProbability().real()
              << ", p|-> = " << sameObj.calculateMinusProbability().real() << std::endl;
    sameObj.evolution(0);
    std::cout << "[T = 0]: " << "p|+> = " << sameObj.calculatePlusProbability().real()
              << ", p|-> = " << sameObj.calculateMinusProbability().real() << std::endl;
    std::cout << "<-------------------------------------------->" << std::endl;

    std::cout << "Measurement: " << sameObj.measure() << std::endl;

    std::cout << "[T = 0]: " << "p|+> = " <<sameObj.calculatePlusProbability().real()
              << ", p|-> = " << sameObj.calculateMinusProbability().real() << std::endl;
    std::cout << "<-------------------------------------------->" << std::endl;

    std::cout << "Evolution for 1 period" << std::endl;
    sameObj.evolution(1);

    std::cout << "[T = 0]: " << "p|+> = " <<sameObj.calculatePlusProbability().real()
              << ", p|-> = " << sameObj.calculateMinusProbability().real() << std::endl;
    std::cout << "<-------------------------------------------->" << std::endl;

    std::cout << "Measurement: " << sameObj.measure() << std::endl;

    std::cout << "[T = 0]: " << "p|+> = " <<sameObj.calculatePlusProbability().real()
              << ", p|-> = " << sameObj.calculateMinusProbability().real() << std::endl;
    std::cout << "<-------------------------------------------->" << std::endl;
}

void testNewObject() {
    std::cout << "Test with NEW object each time" << std::endl;
    std::cout << "First with time 1 period, then 0,5 then 0,25 and finally 0." << std::endl;
    Spin newObj = Spin(1);
    std::cout << "[T = 1]: " << "p|+> = " << newObj.calculatePlusProbability().real()
              << ", p|-> = " << newObj.calculateMinusProbability().real() << std::endl;
    newObj = Spin(0.5);
    std::cout << "[T = 0.5]: " << "p|+> = " << newObj.calculatePlusProbability().real()
              << ", p|-> = " << newObj.calculateMinusProbability().real() << std::endl;
    newObj= Spin(0.25);
    std::cout << "[T = 0.25]: " << "p|+> = " << newObj.calculatePlusProbability().real()
              << ", p|-> = " << newObj.calculateMinusProbability().real() << std::endl;
    newObj= Spin(0);
    std::cout << "[T = 0]: " << "p|+> = " << newObj.calculatePlusProbability().real()
              << ", p|-> = " << newObj.calculateMinusProbability().real() << std::endl;
    std::cout << "<-------------------------------------------->" << std::endl;

    std::cout << "Measurement: " << newObj.measure() << std::endl;

    std::cout << "[T = 0]: " << "p|+> = " << newObj.calculatePlusProbability().real()
              << ", p|-> = " << newObj.calculateMinusProbability().real() << std::endl;
    std::cout << "<-------------------------------------------->" << std::endl;

    std::cout << "New spin at 1 period" << std::endl;
    newObj = Spin(1);

    std::cout << "[T = 0]: " << "p|+> = " << newObj.calculatePlusProbability().real()
              << ", p|-> = " << newObj.calculateMinusProbability().real() << std::endl;
    std::cout << "<-------------------------------------------->" << std::endl;

    std::cout << "Measurement: " << newObj.measure() << std::endl;

    std::cout << "[T = 0]: " << "p|+> = " << newObj.calculatePlusProbability().real()
              << ", p|-> = " << newObj.calculateMinusProbability().real() << std::endl;
    std::cout << "<-------------------------------------------->" << std::endl;
}

void testZeno() {

    std::cout << "Test of Zeno class" << std::endl;
    Zeno* zeno = new Zeno();

    std::cout << "With collapse and 50 repetitions, print turned on." << std::endl;
    zeno->measure(true, true, 0.5, 50);
    std::cout << "<-------------------------------------------->" << std::endl;

    std::cout << "Without collapse and 50 repetitions, print turned on." << std::endl;
    zeno->measure(true, false, 0.5, 50);
    std::cout << "<-------------------------------------------->" << std::endl;

    // measurement at 1 period interval
    zeno->measure(false, true, 1, 100000);
    std::cout << "<-------------------------------------------->" << std::endl;
    // measurement at 0.5 period interval
    zeno->measure(false, true, 0.5, 100000);
    std::cout << "<-------------------------------------------->" << std::endl;
    // measurement at 0.25 period interval;
    zeno->measure(false, true, 0.25, 100000);
    std::cout << "<-------------------------------------------->" << std::endl;

    // simulation of quantum zeno paradox
    zeno->paradox(100000);
}

int main() {

    testSameObject();

    testNewObject();

    testZeno();

    return 0;
}
