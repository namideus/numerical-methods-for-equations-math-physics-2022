import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.theme.Theme;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Objects;

/**
 * @forked by Mirzaali Ruzimatov
 * @author Yiman Altynbek uulu
 * course: Numerical methods for equations of mathematical physics.
 *
 *
 * 2022/11/17
 * */

public class Main extends JFrame{
    // Double
    private static double[] x, array_a, array_b, array_c, array_f, array_r, array_gamma, array_sol, array_sol_origin;
    // Variables
    private static double error, h, t, tau, curant, Tmax, a, b, theta, eta0, eta1, zeta0, zeta1, phi0, phi1, E;
    // Integer
    private static int N = 8, M, T, k, problem = 0, method = 0;
    // X and Y coordinate lists
    private static ArrayList<Double> xData1,yData1,xData2,yData2,xData3,yData3;
    // GUI (Java Swing)
    private static XYChart chart;
    private static XYSeries testFunctionSeries, interpolateFunctionSeries;
    private final JComboBox<Integer> nodesChoice;
    private final JComboBox<String> problemsChoice,methodChoice;
    private final JComboBox<Double> epsilonChoice;
    private final JTextField epsilonInput,kInput;
    private JButton display = new JButton("Display");
    // Series names
    private static final String seriesName1 = "Numerical solution";
    private static final String seriesName2 = "Analytical solution";
    //------------------------------------------------------JFRAME------------------------------------------------------------------
    public Main() {
        setLayout(new BorderLayout());
        JPanel mainPanel = new JPanel();
        add(mainPanel, BorderLayout.CENTER);
        // Parameter, methods selection
        Integer[] choicesNodes = {9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193};
        // Integer[] choices = {9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193};
        Double[] choicesCurant = {0.1,0.5,0.8,1.0,5.0,10.0,20.0,50.0};
        String[] choicesProblem = {"Трансформация «k»-й гармоники ряда Фурье"};
        String[] choicesMethod = {"Явная схема (1)", "Схема Кранка-Николсона (2)", "Схема, сохраняющая монотонность (5)"};
        Double[] choicesEpsilon = { 0.5,0.3,0.1,0.08,0.0625,0.015 };
        JLabel nodesLabel = new JLabel("Nodes");
        nodesChoice = new JComboBox<>(choicesNodes);                    // Node selection
        JLabel problemLabel = new JLabel("Problem");
        problemsChoice = new JComboBox<>(choicesProblem);               // Problems selection
        methodChoice = new JComboBox<>(choicesMethod);                  // Method selection
        JLabel methodLabel = new JLabel("Method");
        epsilonChoice = new JComboBox<>(choicesEpsilon);                // Epsilon selection
        JLabel epsilonLabel = new JLabel("Epsilon");
        epsilonInput = new JTextField("0.0135");
        epsilonInput.setToolTipText("Epsilon");
        JLabel kLabel = new JLabel("k");
        kInput = new JTextField("1");
        kInput.setToolTipText("k");
        //phi1Input = new JTextField("1");
        //JLabel phi1Label = new JLabel("φ1");
        //phi1Input.setToolTipText("Phi 1");
        // Control panel
        JPanel controlPanel = new JPanel();
        controlPanel.setLayout(new GridLayout(40, 1));
        controlPanel.setBackground(Color.lightGray);
        controlPanel.add(nodesLabel);
        controlPanel.add(nodesChoice);
        controlPanel.add(problemLabel);
        controlPanel.add(problemsChoice);
        controlPanel.add(methodLabel);
        controlPanel.add(methodChoice);
        controlPanel.add(epsilonLabel);
        controlPanel.add(epsilonInput);
        controlPanel.add(kLabel);
        controlPanel.add(kInput);
        controlPanel.add(display);
        add(controlPanel, BorderLayout.EAST);

        ActionListener actionListener = e -> {
            try {
                problem = problemsChoice.getSelectedIndex();
                method = methodChoice.getSelectedIndex();
                E = Double.parseDouble(Objects.requireNonNull(epsilonInput.getText()));
                k = Integer.parseInt(Objects.requireNonNull(kInput.getText()));
                N = Integer.parseInt(Objects.requireNonNull(nodesChoice.getSelectedItem()).toString());
                ApplyNumericalMethod();
                Draw(this);
            } catch (InterruptedException ex) {
                throw new RuntimeException(ex);
            }
        };
        problemsChoice.addActionListener(actionListener);
        methodChoice.addActionListener(actionListener);
        nodesChoice.addActionListener(actionListener);
        epsilonInput.addActionListener(actionListener);
        display.addActionListener(actionListener);
    }
    //------------------------------------------------------------------------------------------------------------------------------

    //------------------------------------------------------------------------------------------------------------------------------
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
    public static void addCoord(double x, double y) {
        xData1.add(x);
        yData1.add(y);
    }
    // Chart solutions
    private static void Draw(JFrame jFrame) throws InterruptedException {
        // First clear up array lists
        xData1.clear();
        yData1.clear();
        xData2.clear();
        yData2.clear();
        // Coordinates of numerical solution
        for (int i = 1; i <= N; i++) {
            addCoord(x[i], array_sol[i]);
        }
        // Draw analytical solution
        for(double x=0; x<=1; x+=0.001) {
            xData2.add(x);
            // yData2.add(SolutionFunction(x));
        }
        // Update graphs
        chart.updateXYSeries(seriesName1, xData1, yData1, null);
        chart.updateXYSeries(seriesName2, xData2, yData2, null);
        chart.setTitle("Error: "+error);
        jFrame.repaint();
    }
    //------------------------------------------------------------------------------------------------------------------------------

    //-------------------------------------------------------NUMERICAL-METHODS------------------------------------------------------
    // Double sweep algorithm
    private static double[] DoubleSweepAlgorithm(int n, double[] A, double[] B, double[] C, double[] f) {
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
        for (int i = n-1; i >= 1; i--) {
            U[i] = alpha[i]*U[i+1] + beta[i];
        }
        return U;
    }
    // RHS function
    private static double Function(double x, double t) {
        return 0.0;
    }
    // Analytical solution
    private static double SolutionFunction(double x, double t) {
        return Math.sin(Math.PI*k*x)*Math.pow(Math.E,-(Math.PI*Math.PI*k*k))+x*Psi1(t)+(1-x)*Psi0(t);
    }
    // Phi(x)
    private static double Phi(double x, double t) {
        return Math.sin(Math.PI*k*x)+x*Psi1(t)+(1-x)*Psi0(t);
    }
    // Psi0(t)
    private static double Psi0(double t) {
        return 1.0;
    }
    // Psi1(t)
    private static double Psi1(double t) {
        return 1.0;
    }
    // Error
    public static double Error(double[] a1, double[] a2) {
        double max1=0, max2=0;
        for (int i=1; i<=N; i++) {
            max1 = Math.max(max1, Math.abs(a1[i] - a2[i]));
            max2 = Math.max(max2, Math.abs(a1[i]));
        }
        return max1/max2*100;
    }
    // Numerical method
    private static void ApplyNumericalMethod() {
        // Initial & boundary values
       /* switch (problem) {
            case 0 -> ;
            case 1 -> ;
            case 2 -> ;
        }*/
        // Initialize arrays and variables
        h = 1.0/(N-1);
        tau = Tmax/M;
        curant = (E*tau)/(h*h);

        x = new double[N+1];
        array_a = new double[N+1];
        array_b = new double[N+1];
        array_c = new double[N+1];
        array_f = new double[N+1];
        array_r = new double[N+1];
        array_gamma = new double[N+1];
        array_sol_origin = new double[N+1];

        // Grid
        for (int i = 1; i <= N; i++)
            x[i] = (i-1.0)/(N-1);

        // Analytical solution in the grid
        for (int i = 1; i <= N; i++) {
            //array_sol_origin[i] = SolutionFunction(x[i]);
        }
        // Compute R array
        for (int i = 1; i <= N; i++) {
            //array_r[i] = (CoefficientA(x)*h)/(2*E);
        }
        // Select method
        switch (method) {
            case 0 -> theta = 0;
            case 1 -> theta = 0.5;
            case 2 -> theta = Math.max(0.5, 1-3.0/(4*k));
        }
        // Compute A, B, C, F arrays
        for (int i = 1; i <= N; i++) {
            //x = (i-1.0)/(N-1);
            array_a[i] = E/(h*h)*(array_gamma[i]-array_r[i]);
            //array_b[i] = (2*E*array_gamma[i])/(h*h)+CoefficientB(x);
            array_c[i] = E/(h*h)*(array_gamma[i]+array_r[i]);
           // array_f[i] = -Function(x);
        }
        // Initialize boundary/corner values
        array_a[1] = 0;
        array_b[1] = zeta0+eta0*(E/h);
        array_c[1] = eta0*(E/h);
        array_f[1] = phi0;
        array_a[N] = eta1*(E/h);
        array_b[N] = zeta1+eta1*(E/h);
        array_c[N] = 0;
        array_f[N] = phi1;
        // Call double-sweep algorithm
        array_sol = DoubleSweepAlgorithm(N, array_a, array_b, array_c, array_f);
        // Calculate error
        error = Error(array_sol_origin, array_sol);
    }
    //------------------------------------------------------------------------------------------------------------------------------

    //----------------------------------------------------------MAIN----------------------------------------------------------------
    public static void main(String[] args) throws InterruptedException {
        Setup();
        chart = new XYChartBuilder().width(1750).height(900).xAxisTitle("X").yAxisTitle("Y").build();
        chart.getStyler().setChartBackgroundColor(Color.lightGray);
        chart.getStyler().setCursorBackgroundColor(Color.lightGray);
        chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideS);
        chart.getStyler().setZoomEnabled(true);
        chart.getStyler().setZoomResetByButton(true);
        chart.getStyler().setCursorEnabled(true);
        chart.getStyler().setPlotBorderVisible(true);
        chart.getStyler().setPlotMargin(10);
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
        frame.setTitle("Heat equation");
        frame.add(new XChartPanel<>(chart));
        frame.setSize(frame.getWidth(), frame.getHeight());
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setVisible(true);

        // Testing
        double phase = 0;
        t = 0;
        tau = 2;
        Tmax = 213;
        double[][] initdata = getSineData(phase);
        int timer = 20;
        while (true) {
            phase += 2 * Math.PI * 2 / 20.0;
            Thread.sleep(150);
            final double[][] data = getSineData(phase);
            javax.swing.SwingUtilities.invokeLater(() -> {
                chart.updateXYSeries(seriesName1, data[0], data[1], null);
                frame.repaint();
            });
            t += tau;
            if(t > Tmax) {
                t = 0;
                break;
            }
        }
    }
    private static double[][] getSineData(double phase) {
        double[] xData = new double[100];
        double[] yData = new double[100];
        for (int i = 0; i < xData.length; i++) {
            double radians = phase + (2 * Math.PI / xData.length * i);
            xData[i] = radians;
            yData[i] = Math.sin(radians);
        }
        return new double[][] { xData, yData };
    }
}
//------------------------------------------------------------------------------------------------------------------------------
