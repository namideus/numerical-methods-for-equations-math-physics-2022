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
    private static Main frame;
    private static ComputationThread thread;
    private static double[] x, array_a, array_b, array_c, array_f, array_r, array_gamma, array_sol, array_sol_origin;
    private static double error, h, t, tau, curant, Tmax, a, b, theta, eta0, eta1, zeta0, zeta1, phi0, phi1, E;
    private static int N = 8, M, T, k, problem = 0, scheme = 0;
    private static ArrayList<Double> xData1,yData1,xData2,yData2,xData3,yData3;
    private static final String seriesName1 = "Numerical solution";
    private static final String seriesName2 = "Analytical solution";
    //------------------------------------------------------JFRAME------------------------------------------------------------------
    private static XYChart chart;
    private static XYSeries testFunctionSeries, interpolateFunctionSeries;
    private final JComboBox<Integer> nodesChoice, timeNodesChoice, tMaxChoice, kChoice;
    private final JComboBox<String> problemsChoice,schemeChoice;
    private final JTextField epsilonInput, curantInput, timeInput;
    private JButton computeButton = new JButton("Compute");
    //------------------------------------------------------JFRAME------------------------------------------------------------------
    public Main() {
        setLayout(new BorderLayout());
        JPanel mainPanel = new JPanel();
        add(mainPanel, BorderLayout.CENTER);
        // Parameter, methods selection
        Integer[] choicesNodes = {9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193};
        Integer[] choicesTimeNodes= {5,9,17,33,65,129,257,513};
        Integer[] choicesTmax = {1,2,3,4,5,6,7};
        Integer[] choicesK = {1,2,3,4,5};
        String[] choicesProblem = {"Тестовые задача 1."};
        String[] choicesScheme = {"Явная схема (1)", "Схема Кранка-Николсона (2)", "Схема, сохраняющая монотонность (5)"};
        JLabel nodesLabel = new JLabel("Nodes");
        nodesChoice = new JComboBox<>(choicesNodes);                    // Node selection
        JLabel timeNodesLabel = new JLabel("Time nodes");
        timeNodesChoice = new JComboBox<>(choicesTimeNodes);            // Time node selection
        JLabel tMaxLabel = new JLabel("Tmax");
        tMaxChoice = new JComboBox<>(choicesTmax);                      // Tmax selection
        JLabel kLabel = new JLabel("K");
        kChoice = new JComboBox<>(choicesK);
        JLabel problemLabel = new JLabel("Problem");
        problemsChoice = new JComboBox<>(choicesProblem);               // Problems selection
        JLabel schemeLabel = new JLabel("Scheme");
        schemeChoice = new JComboBox<>(choicesScheme);                  // Method selection
        JLabel epsilonLabel = new JLabel("Epsilon");                // Enter epsilon
        epsilonInput = new JTextField("0.0135");
        JLabel curantLabel = new JLabel("Curant");
        curantInput = new JTextField();
        JLabel timeLabel = new JLabel("Time");
        timeInput = new JTextField();
        // Control panel
        JPanel controlPanel = new JPanel();
        controlPanel.setLayout(new GridLayout(35, 1));
        controlPanel.setBackground(Color.lightGray);
        controlPanel.add(schemeLabel);
        controlPanel.add(schemeChoice);
        controlPanel.add(problemLabel);
        controlPanel.add(problemsChoice);
        controlPanel.add(kLabel);
        controlPanel.add(kChoice);
        controlPanel.add(epsilonLabel);
        controlPanel.add(epsilonInput);
        controlPanel.add(nodesLabel);
        controlPanel.add(nodesChoice);
        controlPanel.add(timeNodesLabel);
        controlPanel.add(timeNodesChoice);
        controlPanel.add(tMaxLabel);
        controlPanel.add(tMaxChoice);
        controlPanel.add(curantLabel);
        controlPanel.add(curantInput);
        controlPanel.add(timeLabel);
        controlPanel.add(timeInput);
        controlPanel.add(computeButton);
        add(controlPanel, BorderLayout.EAST);

        ActionListener actionListener = e -> {
            scheme = schemeChoice.getSelectedIndex();
            problem = problemsChoice.getSelectedIndex();
            k = kChoice.getItemAt(kChoice.getSelectedIndex());
            E = Double.parseDouble(Objects.requireNonNull(epsilonInput.getText()));
            N = nodesChoice.getItemAt(nodesChoice.getSelectedIndex());
            M = timeNodesChoice.getItemAt(timeNodesChoice.getSelectedIndex());
            Tmax = tMaxChoice.getItemAt(tMaxChoice.getSelectedIndex());
            startThread();
        };
        problemsChoice.addActionListener(actionListener);
        schemeChoice.addActionListener(actionListener);
        nodesChoice.addActionListener(actionListener);
        timeNodesChoice.addActionListener(actionListener);
        tMaxChoice.addActionListener(actionListener);
        kChoice.addActionListener(actionListener);
        epsilonInput.addActionListener(actionListener);
        computeButton.addActionListener(actionListener);
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
    private void Draw() throws InterruptedException {
        // Clear array lists
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
       // chart.updateXYSeries(seriesName1, xData1, yData1, null);
        //chart.updateXYSeries(seriesName2, xData2, yData2, null);
        chart.setTitle("Error: "+error);
        repaint();
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
        switch (scheme) {
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
    public static Color color(double val) {
        double H = val * 0.3;
        double S = 0.9;
        double B = 0.9;
        int rgb = Color.HSBtoRGB((float)H, (float)S, (float)B);
        int red = (rgb >> 16) & 0xFF;
        int green = (rgb >> 8) & 0xFF;
        int blue = rgb & 0xFF;
        return new Color(red, green, blue, 0x33);
    }
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
        chart.getStyler().setMarkerSize(5);
        // Series 1
        testFunctionSeries = chart.addSeries(seriesName1, xData1, yData1);
        testFunctionSeries.setLineColor(Color.blue);
        testFunctionSeries.setMarkerColor(color(1.0));
        testFunctionSeries.setLineWidth(1.2f);
        // Series 2
        interpolateFunctionSeries = chart.addSeries(seriesName2, xData2, yData2);
        interpolateFunctionSeries.setLineColor(Color.red);
        interpolateFunctionSeries.setMarkerColor(Color.red);
        interpolateFunctionSeries.setLineWidth(1.2f);
        // Main frame
        frame = new Main();
        frame.setTitle("Heat equation");
        frame.add(new XChartPanel<>(chart));
        frame.setSize(frame.getWidth(), frame.getHeight());
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setVisible(true);
    }
    private void startThread() {
        if(thread == null || thread.isInterrupted())
        {
            thread = new ComputationThread();
        } else if(thread.isAlive())
        {
            thread.interrupt();
            thread = new ComputationThread();
        }
        thread.start();
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

    private class ComputationThread extends Thread {
        public ComputationThread() {
           // new Thread(this).start();
        }
        @Override
        public void run() {
            try {
                double phase = 0;
                t = 0;
                tau = 2;
                Tmax = 100;
                while (t <= Tmax) {
                    phase += 2 * Math.PI * 2 / 20.0;
                    sleep(130);
                    double[][] data = getSineData(phase);
                    chart.updateXYSeries(seriesName1, data[0], data[1], null);
                    frame.repaint();
                    t += tau;

                  /*  ApplyNumericalMethod();
                    Draw();*/
                }
                interrupt();
            } catch(Exception e) {
                System.out.println("Exception is caught: " + e);
            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------------
