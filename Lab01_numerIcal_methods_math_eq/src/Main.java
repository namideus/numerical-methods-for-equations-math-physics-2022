import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Objects;

/**
 * @author Yiman Altynbek uulu
 *
 * course: Numerical methods for equations of mathematical physics.
 *
 *
 * 2022/09/28
 * */

public class Main extends JFrame{
    // Double
    private static double[] arr, arr_sol;
    // Variables
    private static double u0, T = 1,h, a = 0.5, lambda = 0.5, error;
    private static double theta, eta0 = 0, eta1 = 1, zeta0 = 0, zeta1 = 1, phi0 = 0, phi1 = 1, E;
    // Integer
    private static int N = 8, problem = 0, method = 0; // M:1,2;
    // X and Y coordinate lists
    private static ArrayList<Double> xData1;
    private static ArrayList<Double> yData1;
    private static ArrayList<Double> xData2;
    private static ArrayList<Double> yData2;
    private static ArrayList<Double> xData3;
    private static ArrayList<Double> yData3;
    // User Interface (Java Swing)
    private static XYChart chart;
    private static XYSeries testFunctionSeries, interpolateFunctionSeries;
    private final JComboBox<Integer> nodesChoice;
    private final JComboBox<String> problemsChoice;
    private final JComboBox<String> methodChoice;
    private final JComboBox<Double> epsilonChoice;

    private final JTextField epsilonInput;
    private final JTextField thetaInput;
    private final JTextField etaInput;

    private JButton display = new JButton("Display");
    // Series names
    private static final String seriesName1 = "Test function";
    private static final String seriesName2 = "Analytical solution";
    private static final String seriesName3 = "Whatever";

    // Building the user interface
    public Main() {
        setLayout(new BorderLayout());
        JPanel mainPanel = new JPanel();
        mainPanel.setBackground(Color.WHITE);
        add(mainPanel, BorderLayout.CENTER);
        // Parameter, methods selection
        // Customize parameters!!
        Integer[] choicesNodes = { 8, 16, 32, 64, 128, 256, 512, 1024 };
        String[] choicesProblem = { "Problem 1", "Problem 5", "Problem 13.1", "Problem 13.4" };
        String[] choicesMethod = { "Method 1", "Method 2" };
        Double[] choicesEpsilon = { 0.5, 0.3, 0.1, 0.08, 0.0625, 0.015 };

        nodesChoice = new JComboBox<>(choicesNodes);                    // Node selection
        problemsChoice = new JComboBox<>(choicesProblem);               // Problems selection
        methodChoice = new JComboBox<>(choicesMethod);                  // Method selection
        epsilonChoice = new JComboBox<>(choicesEpsilon);                // Epsilon selection

        epsilonInput = new JTextField();                                // Epsilon input
        thetaInput = new JTextField();
        etaInput = new JTextField();
        // Control panel
        JPanel controlPanel = new JPanel();
        controlPanel.setLayout(new GridLayout(8, 1));
        controlPanel.setBackground(Color.LIGHT_GRAY);
        controlPanel.add(nodesChoice);
        controlPanel.add(problemsChoice);
        controlPanel.add(methodChoice);
        controlPanel.add(epsilonChoice);
        controlPanel.add(epsilonInput);
        controlPanel.add(thetaInput);
        controlPanel.add(etaInput);
        controlPanel.add(display);
        add(controlPanel, BorderLayout.EAST);

        // Listen to display button click, update the graph
        display.addActionListener(actionEvent -> {
            N = Integer.parseInt(Objects.requireNonNull(nodesChoice.getSelectedItem()).toString());             // Get number of nodes
            problem = problemsChoice.getSelectedIndex();                                                        // Get problem ubdex
            method = methodChoice.getSelectedIndex();                                                           // Get method index
            E = Double.parseDouble(Objects.requireNonNull(epsilonChoice.getSelectedItem()).toString());         // Get epsilon index
            Draw();                                                                                             // Draw graphs
            repaint();                                                                                          // Show the change
        });
    }

    /**
     * Algorithm of "progonka"
     * */
    private double[] Progonka(int n, double[] A, double[] B, double[] C, double[] f)
    {
        double[] U = new double[n+1];
        double[] alpha = new double[n+1];
        double[] beta = new double[n+1];

        alpha[1] = C[1]/B[1];
        beta[1] = f[1]/B[1];

        for (int i = 2; i <= n-1; i++) {
            alpha[i] = C[i]/(B[i] - alpha[i-1]*A[i]);
            beta[i] = (f[i] + beta[i-1]*A[i])/(B[i] - alpha[i-1]*A[i]);
        }

        U[n] = (f[n] + beta[n-1]*A[n])/(B[n] - alpha[n-1]*A[n]);

        for (int i = n-1; i >= 1; i--)
            U[i] = alpha[i]*U[i+1] + beta[i];

        return U;
    }
    /**
     * Test problems
     *
     * P: 1, 5, 13.1, 13.4
     *
     * */
    private static double Function(double x, double y) {
        // Change this function
        return ((1 - 2.0*x)/(1.0 + x) - E - Math.pow(Math.E, -1.0 / E)/(1.0 - Math.pow(Math.E, -1.0 / E)))
                * (2.0 / Math.pow(1.0 + x, 3));
    }
    /**
     * Analytical solutions.
     *
     * P: 1, 2, 4
     *
     * */
    private static double SolutionFunction(double x) {
        return switch (problem) {
            case 0 -> (phi0/zeta0 - (phi1/zeta1 - phi0/zeta0 + E - 0.5)/(Math.pow(Math.E, -1/E)-1))
                    +Math.pow(Math.E, -x/E)*((phi1/zeta0 - phi0/zeta0 + E - 0.5)/(Math.pow(Math.E, -1/E)-1))
                    -E*x + 0.5 * x * x;
            case 1 -> (phi1/zeta1 + E - 0.5 - Math.pow(Math.E, -1/E)*(phi0-E*E)/(1+zeta0))*
                    (1+zeta0)/(1+zeta0-Math.pow(Math.E, -1/E)*zeta0) +
                    Math.pow(Math.E, -x/E)*((phi0-E*E)/(1+zeta0)-(zeta0)/(1+zeta0)*(phi1/zeta1
                    +E-0.5-Math.pow(Math.E, -1/E)*(phi0-E*E)/(1+zeta0))*(1/(1+zeta0-Math.pow(Math.E, -1/E)*zeta0)))
                    -E*x + 0.5 * x * x;
            case 2 -> x / (1.0 + x) + (Math.pow(Math.E, -1.0 / E) - Math.pow(Math.E,
                    -(2.0 * x) / (E * (1.0 + x)))) / (2.0 * (1 - Math.pow(Math.E, -1.0 / E)));
            default -> 0;
        };
    }

    /**
     * Compute error
     * */
    public static double err(double[] m1, double[] m2)
    {
        // Change this function
        return 0;
    }

    // Set up
    private static void Setup() {
        xData1 = new ArrayList<>();
        yData1 = new ArrayList<>();
        xData2 = new ArrayList<>();
        yData2 = new ArrayList<>();
        xData1.add(0d);
        yData1.add(0d);
        xData2.add(0d);
        yData2.add(0d);
    }

    // Draw
    private static void Draw() {
        // First clear up array lists
        xData1.clear();
        yData1.clear();
        xData2.clear();
        yData2.clear();
        //---------------------------------------------------------------
        // Initial values set up (or remove it)

     /*   h = T/(N-1);
        switch(problem) {
            case 0:
                u0 = 2; break;
            case 1:
                u0 = 0; break;
            case 2:
                u0 = 1; break;
            case 3:
                u0 = phi;
        }*/

        // Draw graph of solution function
        //        for (double x = 0.0; x <= 1.0; x += 0.001) {
        //            xData1.add(x);
        //            yData1.add(SolutionFunc(x));
        //        }
        // Initialize certain variables here ...

        // Solve by a selected method
        switch(method) {
            // Method 1
            case 0:
                for(double x=0; x<=10; x+=0.01) {
                    xData2.add(x);
                    yData2.add(SolutionFunction(x)); // Just a sample
                }
                break;
            // Method 2
            case 1:
                for(double x=0; x<=10; x+=0.01) {
                    xData1.add(x);
                    yData1.add(Math.cos(x)); // Just a sample
                }
                break;
        }
        // Update graphs
        chart.updateXYSeries(seriesName1, xData1, yData1, null);
        chart.updateXYSeries(seriesName2, xData2, yData2, null);

        // Compute the error and display, redefine error function
        error = err(arr, arr_sol);
        chart.setTitle("error: "+error);
    }

    // Main
    public static void main(String[] args) {
        Setup();
        // Create Chart
        chart = new XYChartBuilder().width(1800).height(1000).xAxisTitle("X").yAxisTitle("Y").build();
        // Customize Chart
        chart.getStyler().setChartBackgroundColor(Color.LIGHT_GRAY);
        chart.getStyler().setCursorBackgroundColor(Color.GRAY);
        chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideS);
        chart.getStyler().setZoomEnabled(true);
        chart.getStyler().setZoomResetByButton(true);
        chart.getStyler().setMarkerSize(0);

        // Series 1
        testFunctionSeries = chart.addSeries(seriesName1, xData1, yData1);
        testFunctionSeries.setLineColor(Color.blue);
        testFunctionSeries.setLineWidth(1.2f);

        // Series 2
        interpolateFunctionSeries = chart.addSeries(seriesName2, xData2, yData2);
        interpolateFunctionSeries.setLineColor(Color.RED);
        interpolateFunctionSeries.setLineWidth(1.2f);

        Main frame = new Main();
        frame.setTitle("Project #1");
        frame.add(new XChartPanel<>(chart));
        frame.setSize(frame.getWidth(), frame.getHeight());
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setVisible(true);
    }
}
