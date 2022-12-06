import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.theme.MatlabTheme;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Objects;

/**
 * @author Yiman Altynbek uulu
 * course: Numerical methods for equations of mathematical physics.
 *
 *
 * 2022/11/17
 * */

public class Main extends JFrame{
    private static Main frame;
    private static ComputationThread thread;
    private static double[] x, alpha, beta, array_u, array_v, array_a, array_b, array_c, array_f, array_sol_origin, U;
    private static double error=0.0, h, t=0, tau, curant, Tmax, theta, E, A, B, C, A0, B0, C0;
    private static boolean firstCycle = true, lastCycle = true;
    private static int N = 8, M, m = 1, T, k, problem = 0, scheme = 0;
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
    private JButton exitButton = new JButton("Exit");
    //------------------------------------------------------JFRAME------------------------------------------------------------------
    public Main() {
        setLayout(new BorderLayout());
        JPanel mainPanel = new JPanel();
        add(mainPanel, BorderLayout.CENTER);
        // Parameter, methods selection
        Integer[] choicesNodes = {5,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193};
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
        controlPanel.setBackground(Color.white);
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
        controlPanel.add(exitButton);
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
        exitButton.addActionListener(e -> System.exit(0));
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
    // Graph solutions
    private void Graph() {
        // Clear array lists
        xData1.clear();
        yData1.clear();
        xData2.clear();
        yData2.clear();
        // Numerical solution coordinates
        for (int i = 1; i <= N; i++) {
            xData1.add(x[i]);
            yData1.add(array_v[i]);
        }
        // Analytical solution coordinates
        for(int i = 1; i <= 100; i++) {
            double xe = (i - 1.0) / (100 - 1.0);
            xData2.add(xe);
            yData2.add(U(xe, t));
        }
        // Update graphs
        chart.updateXYSeries(seriesName1, xData1, yData1, null);
        chart.updateXYSeries(seriesName2, xData2, yData2, null);
        chart.setTitle("Error: "+error);
        repaint();
    }
    //------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------NUMERICAL-METHODS------------------------------------------------------
    // Double sweep algorithm
    private static double[] ConstantDoubleSweep(int N, double A, double B, double C, double[] F) {
        alpha = new double[N+1];
        beta = new double[N+1];
        U = new double[N+1];
        alpha[1] = Math.floor(C/B);
        beta[1] = F[1]/B;
        for (int i = 2; i <= N-1; i++) {
            double v = B - alpha[i - 1] * A;
            alpha[i] = C / v;
            beta[i] = (F[i] + beta[i - 1]*A) / v;
        }
        U[N] = Math.ceil((F[N] + beta[N-1]*A)/(B - alpha[N-1]*A));
        for (int i = N-1; i >= 1; i--) {
            U[i] = alpha[i]*U[i+1] + beta[i];
        }
        return U;
    }
    // Function
    private static double Function(double x, double t) {
        return 0.0;
    }
    // Analytical solution
    private static double U(double x, double t) {
        return Math.sin(Math.PI*k*x)*Math.pow(Math.E,-(Math.PI*Math.PI*k*k)*E*t)+x*Psi1(t)+(1-x)*Psi0(t);
    }
    // Phi(x)
    private static double Phi(double x, double t) {
        return Math.sin(Math.PI*k*x)+x*Psi1(t)+(1-x)*Psi0(t);
    }
    // Psi0(t)
    private static double Psi0(double t) {
        return 0;
    }
    // Psi1(t)
    private static double Psi1(double t) {
        return 1;
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
    private void ApplyNumericalMethod() {
        // Initialization
        h = 1.0/(N-1);
        tau = Tmax/M;
        curant = (E*tau)/(h*h);
        curantInput.setText(Double.toString(curant));
        timeInput.setText(Double.toString(t));
        // Select a scheme
        switch (scheme) {
            case 0 -> theta = 0;
            case 1 -> theta = 1/2.0;
            case 2 -> theta = Math.max(0.5, 1-3/(4.0*curant));
        }
        // Constants
        A = theta*curant;
        C = theta*curant;
        B = A + C + 1;
        A0 = (1-theta)*curant;
        C0 = (1-theta)*curant;
        B0 = 1 - A0 - C0;
        // Arrays
        x = new double[N+1];
        array_f = new double[N+1];
        array_u = new double[N+1];
        array_sol_origin = new double[N+1];
        // Grid
        for (int i = 1; i <= N; i++) {
            x[i] = (i-1.0)/(N-1);
        }
        // (V.2.8)-(V.2.10)
        if(t==0) {
            for (int i = 1; i <=N; i++) {
                array_f[i] = Phi(x[i], t);
            }
            array_v = array_f;
        } else {
            // Compute A, B, C, F arrays
            array_f[1] = Psi0(t);
            array_f[N] = Psi1(t);

            for (int i = 2; i < N; i++) {
                array_f[i]=A0*array_v[i-1]+B0*array_v[i]+C0*array_v[i+1];
            }

            array_v = ConstantDoubleSweep(N, A, B, C, array_f);
        }
        // Calculate error
        for (int i = 1; i <= N; i++) {
            array_sol_origin[i] = U(x[i], t);
        }
        error = Math.max(Error(array_sol_origin, array_v), error);
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
        chart.getStyler().setChartBackgroundColor(Color.white);
        chart.getStyler().setTheme(new MatlabTheme());
        chart.getStyler().setCursorBackgroundColor(Color.lightGray);
        chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideNW);
        chart.getStyler().setZoomEnabled(true);
        chart.getStyler().setZoomResetByButton(true);
        chart.getStyler().setCursorEnabled(true);
        chart.getStyler().setPlotBorderVisible(true);
        chart.getStyler().setPlotMargin(10);
        chart.getStyler().setMarkerSize(4);
        // Series 1
        testFunctionSeries = chart.addSeries(seriesName1, xData1, yData1);
        testFunctionSeries.setLineColor(Color.blue);
        testFunctionSeries.setMarkerColor(Color.blue);
        testFunctionSeries.setLineWidth(1.2f);
        // Series 2
        interpolateFunctionSeries = chart.addSeries(seriesName2, xData2, yData2);
        interpolateFunctionSeries.setLineColor(Color.red);
        interpolateFunctionSeries.setMarkerColor(color(1.0));
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
    private class ComputationThread extends Thread {
        public ComputationThread() {}
        @Override
        public void run() {
            try {
                t = 0;
                m = 0;
                error = 0;
                while (t <= Tmax) {
                    sleep(170);
                    ApplyNumericalMethod();
                    Graph();
                    t += tau;
                }
                interrupt();
            } catch(Exception e) {
                System.out.println("Exception is caught: " + e);
            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------------