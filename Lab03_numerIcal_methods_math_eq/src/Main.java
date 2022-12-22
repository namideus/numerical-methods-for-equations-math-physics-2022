import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.lines.MatlabSeriesLines;
import org.knowm.xchart.style.theme.GGPlot2Theme;
import org.knowm.xchart.style.theme.MatlabTheme;
import org.knowm.xchart.style.theme.Theme;
import org.knowm.xchart.style.theme.XChartTheme;

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
    private static double[] x, y;
    private static double[][] u, v, uprev;
    private static double l,m, error=0.0, h, tau, theta, sigma;
    private static int N = 8, method = 0, iteration = 0, problem, k;
    private static ArrayList<Double> xData1,yData1,xData2,yData2,xData3,yData3,xData4,yData4;
    private static final String numericalSeriesName = "Numerical solution";
    private static final String analyticalSeriesName = "Analytical solution";
    //------------------------------------------------------JFRAME------------------------------------------------------------------
    private static XYChart chart1, chart2;
    private static XYSeries solutionFunctionSeries1, solutionFunctionSeries2, numericalFunctionSeries1, numericalFunctionSeries2;
    private final JComboBox<Integer> nodesChoice;
    private final JComboBox<String> problemsChoice, iterationChoice;
    private final JTextField sigmaInput, iterationsInput, timeInput, errorInput;
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
        String[] choicesProblem = {"Problem No: 1 (l=2, m=1)"};
        String[] choicesIterationMethod = {"Jacob", "Seidel"};
        JLabel nodesLabel = new JLabel("Nodes");
        nodesChoice = new JComboBox<>(choicesNodes);                    // Node selection
        JLabel problemLabel = new JLabel("Problem");
        problemsChoice = new JComboBox<>(choicesProblem);               // Problems selection
        JLabel iterationLabel = new JLabel("Iteration method");
        iterationChoice = new JComboBox<>(choicesIterationMethod);      // Iteration selection
        JLabel iterationsLabel = new JLabel("Iterations (K)");
        iterationsInput = new JTextField("1");
        JLabel sigmaLabel = new JLabel("Sigma");                   // Enter sigma
        sigmaInput = new JTextField("0.0135");
        JLabel errorLabel = new JLabel("Error");
        errorInput = new JTextField();
        JLabel timeLabel = new JLabel("Time");
        timeInput = new JTextField();
        // Control panel
        JPanel controlPanel = new JPanel();
        controlPanel.setLayout(new GridLayout(35, 1));
        controlPanel.setBackground(Color.white);
        controlPanel.add(iterationLabel);
        controlPanel.add(iterationChoice);
        controlPanel.add(problemLabel);
        controlPanel.add(problemsChoice);
        controlPanel.add(nodesLabel);
        controlPanel.add(nodesChoice);
        controlPanel.add(iterationsLabel);
        controlPanel.add(iterationsInput);
        controlPanel.add(sigmaLabel);
        controlPanel.add(sigmaInput);
        controlPanel.add(errorLabel);
        controlPanel.add(errorInput);
        controlPanel.add(computeButton);
        controlPanel.add(exitButton);
        add(controlPanel, BorderLayout.EAST);

        ActionListener actionListener = e -> {
            method = iterationChoice.getSelectedIndex();
            problem = problemsChoice.getSelectedIndex();
            N = nodesChoice.getItemAt(nodesChoice.getSelectedIndex());
            k = Integer.parseInt(Objects.requireNonNull(iterationsInput.getText()));
            sigma = Double.parseDouble(Objects.requireNonNull(sigmaInput.getText()));
            startThread();
        };
        problemsChoice.addActionListener(actionListener);
        iterationChoice.addActionListener(actionListener);
        nodesChoice.addActionListener(actionListener);
        iterationsInput.addActionListener(actionListener);
        sigmaInput.addActionListener(actionListener);
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
        xData3 = new ArrayList<>();
        yData3 = new ArrayList<>();
        xData4 = new ArrayList<>();
        yData4 = new ArrayList<>();
        xData1.add(0d);
        yData1.add(0d);
        xData2.add(0d);
        yData2.add(0d);
        xData3.add(0d);
        yData3.add(0d);
        xData4.add(0d);
        yData4.add(0d);
    }
    // Graph solutions
    private void Graph() {
        // Clear array lists
        xData1.clear();
        yData1.clear();
        xData2.clear();
        yData2.clear();
        xData3.clear();
        yData3.clear();
        xData4.clear();
        yData4.clear();
        // Numerical solution
        for (int i = 1; i <= N; i++) {
            xData1.add(x[i]);
            yData1.add(u[i][(int)Math.ceil(N/2.0)]);
        }
        for (int j = 1; j <= N; j++) {
            xData3.add(y[j]);
            yData3.add(u[(int)Math.ceil(N/4.0)][j]);
        }
        // Analytical solution
        for(int i = 1; i <= 100; i++) {
            double xe = (i - 1.0) / (100 - 1.0);
            xData2.add(xe);
            yData2.add(U(xe, 0.5));
        }
        for(int j = 1; j <= 100; j++) {
            double ye = (j - 1.0) / (100 - 1.0);
            xData4.add(ye);
            yData4.add(U(0.25, ye));
        }
        // Update charts
        chart1.updateXYSeries(numericalSeriesName, xData1, yData1, null);
        chart1.updateXYSeries(analyticalSeriesName, xData2, yData2, null);
        chart2.updateXYSeries(numericalSeriesName, xData3, yData3, null);
        chart2.updateXYSeries(analyticalSeriesName, xData4, yData4, null);
        repaint();
    }
    //------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------NUMERICAL-METHODS------------------------------------------------------
    // Function
    private static double F(int i, int j) {
        double lambda=4.0/(h*h)*(Math.sin(Math.PI*l*h/2)*Math.sin(Math.PI*l*h/2)+Math.sin(Math.PI*m*h/2)*Math.sin(Math.PI*m*h/2));
        return lambda*Math.sin(Math.PI*l*x[i])*Math.sin(Math.PI*m*y[j]);
    }
    // Analytical solution
    private static double U(double x, double y) {
        return Math.sin(Math.PI*l*x)*Math.sin(Math.PI*m*y);
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
    // Initialization
    private void Initialize() {
        l = 2;
        m = 1;
        h = 1.0/(N-1);
       // sigma = 1.0/Math.pow(N-1, 3);
        sigmaInput.setText(Double.toString(sigma));
        switch (method) {
            case 0 -> {
                theta = 0;
                tau = 1;
            }
            case 1 -> {
                theta = 1.0;
                tau = 1.0;
            }
        }
        // Arrays
        x = new double[N+1];
        y = new double[N+1];
        v = new double[N+1][N+1];
        u = new double[N+1][N+1];
        uprev = new double[N+1][N+1];

        // Grid
        for (int i = 1; i <= N; i++) {
            x[i] = (i-1.0)*h;
        }
        for (int j = 1; j <= N; j++) {
            y[j] = (j-1.0)*h;
        }
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                u[i][j] = v[i][j] = 0.0;
            }
        }
    }
    // Numerical method
    private void ApplyNumericalMethod() {
        for (int j = 2; j < N; j++) {
            for (int i = 2; i < N; i++) {
                u[i][j] = (theta/4.0)*(u[i-1][j]+u[i][j-1])
                        +(tau-theta)/4.0*(v[i-1][j]+v[i][j-1])
                        +tau/4.0*(v[i+1][j]+v[i][j+1])
                        +(1-tau)*v[i][j]
                        +(tau*h*h)/4.0*F(i, j);
            }
        }
        v = u; // Copy previous numerical solution
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
    public static void main(String[] args) {
        Setup();
        chart1 = new XYChartBuilder().width(850).height(850).xAxisTitle("X").yAxisTitle("Y").build();
        chart1.getStyler().setChartBackgroundColor(Color.lightGray);
        chart1.getStyler().setTheme(new MatlabTheme());
        chart1.getStyler().setLegendPosition(Styler.LegendPosition.InsideNW);
        chart1.getStyler().setZoomEnabled(true);
        chart1.getStyler().setZoomResetByButton(true);
        chart1.getStyler().setPlotBorderVisible(true);
        chart1.getStyler().setPlotMargin(10);
        chart1.getStyler().setMarkerSize(4);
        chart2 = new XYChartBuilder().width(850).height(850).xAxisTitle("X").yAxisTitle("Y").build();
        chart2.getStyler().setChartBackgroundColor(Color.lightGray);
        chart2.getStyler().setTheme(new MatlabTheme());
        chart2.getStyler().setLegendPosition(Styler.LegendPosition.InsideNW);
        chart2.getStyler().setZoomEnabled(true);
        chart2.getStyler().setZoomResetByButton(true);
        chart2.getStyler().setPlotBorderVisible(true);
        chart2.getStyler().setPlotMargin(10);
        chart2.getStyler().setMarkerSize(4);
        // Chart 1
        solutionFunctionSeries1 = chart1.addSeries(numericalSeriesName, xData1, yData1);
        solutionFunctionSeries1.setLineColor(Color.blue);
        solutionFunctionSeries1.setMarkerColor(Color.blue);
        solutionFunctionSeries1.setLineWidth(1.2f);

        numericalFunctionSeries1 = chart1.addSeries(analyticalSeriesName, xData2, yData2);
        numericalFunctionSeries1.setLineColor(Color.red);
        numericalFunctionSeries1.setMarkerColor(color(1.0));
        numericalFunctionSeries1.setLineWidth(1.2f);
        // Chart 2
        solutionFunctionSeries2 = chart2.addSeries(numericalSeriesName, xData3, yData3);
        solutionFunctionSeries2.setLineColor(Color.blue);
        solutionFunctionSeries2.setMarkerColor(Color.blue);
        solutionFunctionSeries2.setLineWidth(1.2f);

        numericalFunctionSeries2 = chart2.addSeries(analyticalSeriesName, xData4, yData4);
        numericalFunctionSeries2.setLineColor(Color.red);
        numericalFunctionSeries2.setMarkerColor(color(1.0));
        numericalFunctionSeries2.setLineWidth(1.2f);
        // Main frame
        frame = new Main();
        frame.setLayout(new FlowLayout());
        frame.setTitle("Poisson equation");
        frame.getContentPane().add(new XChartPanel<>(chart1));
        frame.getContentPane().add(new XChartPanel<>(chart2));
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
    private double chebyshevNormDifference(double[][] u)
    {
        double res = 0, sum, max = 0;
        for (int j = 1; j <= N; j++) {
            sum = 0;
            for (int i = 1; i <= N; i++) {
                // max = Math.max(max, Math.abs(u[i][j]-uprev[i][j]));
                sum += Math.abs(Math.abs(u[j][i])-Math.abs(uprev[j][i]));
            }
            res = Math.max(sum, res);
        }
        uprev = u;
        return res;
    }
    private class ComputationThread extends Thread {
        public ComputationThread() {}
        @Override
        public void run() {
            try {
                iteration = 0;
                error = 0;
                errorInput.setText("");
                Initialize();
               //while (iteration <= k) {
                   do {
                       // iterationsInput.setText(iteration+"");
                       sleep(100);
                       ApplyNumericalMethod();
                       Graph();
                       iteration++;
                   }
               //}
               while(chebyshevNormDifference(u) > sigma);
                // Calculate error
                for (int j = 1; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        error = Math.max(Math.abs(u[i][j]-U(x[i], y[j])), error);
                    }
                }
                errorInput.setText(error+"");
                interrupt();
            } catch(Exception e) {
                System.out.println("Exception is caught: " + e);
            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------------