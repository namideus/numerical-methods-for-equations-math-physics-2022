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
    private static double[] array_a, array_b, array_f, array_r, array_gamma, array_sol;
    // Variables
    private static double u0, T = 1.0, h, error, x;
    private static double a, b, theta, eta0, eta1, zeta0, zeta1, phi0, phi1, E;
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
    private final JTextField zeta0Input;
    private final JTextField zeta1Input;
    private final JTextField eta0Input;
    private final JTextField eta1Input;
    private JButton display = new JButton("Display");
    // Series names
    private static final String seriesName1 = "Numerical solution";
    private static final String seriesName2 = "Analytical solution";

    // Building the user interface
    public Main() {
        setLayout(new BorderLayout());
        JPanel mainPanel = new JPanel();
        mainPanel.setBackground(Color.WHITE);
        add(mainPanel, BorderLayout.CENTER);
        // Parameter, methods selection
        // Customize parameters!!
        Integer[] choicesNodes = { 8, 16, 32, 64, 128, 256, 512, 1024 };
        String[] choicesProblem = { "Problem 1", "Problem 2", "Problem 4" };
        String[] choicesMethod = { "With central difference", "With directional difference" };
        Double[] choicesEpsilon = { 0.5, 0.3, 0.1, 0.08, 0.0625, 0.015 };

        nodesChoice = new JComboBox<>(choicesNodes);                    // Node selection
        problemsChoice = new JComboBox<>(choicesProblem);               // Problems selection
        methodChoice = new JComboBox<>(choicesMethod);                  // Method selection
        epsilonChoice = new JComboBox<>(choicesEpsilon);                // Epsilon selection

        epsilonInput = new JTextField();
        epsilonInput.setToolTipText("Epsilon");
        thetaInput = new JTextField();
        thetaInput.setToolTipText("Theta");
        zeta0Input = new JTextField();
        zeta0Input.setToolTipText("Zeta 0");
        zeta1Input = new JTextField();
        zeta1Input.setToolTipText("Zeta 1");
        eta0Input = new JTextField();
        eta0Input.setToolTipText("Eta 0");
        eta1Input = new JTextField();
        eta1Input.setToolTipText("Eta 1");

        // Control panel
        JPanel controlPanel = new JPanel();
        controlPanel.setLayout(new GridLayout(25, 1));
        controlPanel.setBackground(Color.LIGHT_GRAY);
        controlPanel.add(nodesChoice);
        controlPanel.add(problemsChoice);
        controlPanel.add(methodChoice);
        controlPanel.add(epsilonInput);
      /*  controlPanel.add(zeta0Input);
        controlPanel.add(zeta1Input);
        controlPanel.add(eta0Input);
        controlPanel.add(eta1Input);*/
        controlPanel.add(thetaInput);
        controlPanel.add(display);
        add(controlPanel, BorderLayout.EAST);

        // Listen to display button click, update the graph
        display.addActionListener(actionEvent -> {
            E = Double.parseDouble(Objects.requireNonNull(epsilonInput.getText()));
           /* zeta0 = Double.parseDouble(Objects.requireNonNull(zeta0Input.getText()));
            zeta1 = Double.parseDouble(Objects.requireNonNull(zeta1Input.getText()));
            eta0 = Double.parseDouble(Objects.requireNonNull(eta0Input.getText()));
            eta1 = Double.parseDouble(Objects.requireNonNull(eta1Input.getText()));*/
            theta = Double.parseDouble(Objects.requireNonNull(thetaInput.getText()));
            N = Integer.parseInt(Objects.requireNonNull(nodesChoice.getSelectedItem()).toString());             // Get number of nodes
            problem = problemsChoice.getSelectedIndex();                                                        // Get problem ubdex
            method = methodChoice.getSelectedIndex();                                                           // Get method index
            Draw();
            repaint();
        });
    }

    // Algorithm of "progonka"
    private static double[] Progonka(int n, double[] A, double[] B, double[] C, double[] f) {
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

    // Test problems
    private static double Function(double x) {
        return switch (problem) {
            case 0, 1 -> x;
            case 2 -> ((1-2*x)/(1+x)-E-Math.pow(Math.E,-1/E)/(1-Math.pow(Math.E,-1/E)))*(2/Math.pow(1+x,3));
            default -> 0;
        };
    }

    // Analytical solutions.
    private static double SolutionFunction(double x) {
        return switch (problem) {
            case 0 -> phi0+((phi1-phi0+E-0.5)*(Math.pow(Math.E, -x/E)-1))/(Math.pow(Math.E, -1/E)-1)-E*x+0.5*x*x;
            case 1 -> phi0-E*E+((phi0-phi1-E-E*E+0.5)*(Math.pow(Math.E,-x/E)-2))/(2-Math.pow(Math.E,-1/E)-E*x+0.5*x*x);
            case 2 -> x/(1+x)+(Math.pow(Math.E,-1/E)-Math.pow(Math.E,-(2*x)/(E*(1+x))))/(2*(1-Math.pow(Math.E,-1/E)));
            default -> 0;
        };
    }

    // Coefficient "a"
    private static double CoefficientA(double x) {
        return switch (problem) {
            case 0, 1 -> 1.0;
            case 2 -> 2.0/((1+x)*(1+x));
            default -> 0;
        };
    }

    // Coefficient "b"
    private static double CoefficientB(double x) {
        return switch (problem) {
            case 0, 1 -> 0.0;
            case 2 -> 4.0/((1+x)*(1+x)*(1+x));
            default -> 0;
        };
    }

    // Compute error
    public static double Error(double[] m1, double[] m2) {
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
        // Initial & boundary values
        switch (problem) {
            case 0 -> {
                zeta0 = zeta1 = 1.0;
                eta0 = eta1 = 0.0;
                phi0 = 0;
                phi1 = 1;
            }
            case 1 -> {
                zeta0 = zeta1 = 1.0;
                eta0 = 1.0;
                eta1 = 0.0;
                phi0 = 0;
                phi1 = 1;
            }
            case 2 -> {
                zeta0 = zeta1 = 2.0;
                eta0 = 1.0;
                eta1 = 4.0;
                phi0 = -(1+E+1/(1-Math.pow(Math.E,-1/E)));
                phi1 = 1+E+Math.pow(Math.E,-1/E)/(1-Math.pow(Math.E,-1/E));
            }
        }

        h = T/(N-1);
        array_a = new double[N+1];
        array_b = new double[N+1];
        array_f = new double[N+1];
        array_r = new double[N+1];

        for (int i = 1; i <= N; i++) {
            x = (i-1.0)/(N-1);
            array_a[i] = CoefficientA(x);
            array_b[i] = CoefficientB(x);
            array_f[i] = Function(x);
            array_r[i] = (array_a[i]*h)/(2*E);
        }

        // Draw analytical solution
        for(double x=0; x<=1; x+=0.001) {
            xData2.add(x);
            yData2.add(SolutionFunction(x));
        }

        //array_sol = Progonka(N, array_a, array_b, array_f)

        switch(method) {
            // Method 1
            case 0:
                for (int i = 1; i <= N; i++) {
                    array_gamma[i] = 1.0;
                }
                break;
            // Method 2
            case 1:
                for (int i = 1; i <= N; i++) {
                    array_gamma[i] = 1 + Math.abs(array_r[i]);
                }
                break;
        }
        // Update graphs
        chart.updateXYSeries(seriesName1, xData1, yData1, null);
        chart.updateXYSeries(seriesName2, xData2, yData2, null);

        // Compute the error and display, redefine error function
        // error = Error(arr, arr_sol);
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
        // Main frame
        Main frame = new Main();
        frame.setTitle("Advection-Diffusion");
        frame.add(new XChartPanel<>(chart));
        frame.setSize(frame.getWidth(), frame.getHeight());
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setVisible(true);
    }
}
