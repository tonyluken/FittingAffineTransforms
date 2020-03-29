import java.awt.Color;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.BasicStroke; 

import org.jfree.chart.ChartPanel; 
import org.jfree.chart.JFreeChart; 
import org.jfree.data.xy.XYDataset; 
import org.jfree.data.xy.XYSeries; 
import org.jfree.ui.ApplicationFrame; 
import org.jfree.ui.RefineryUtilities; 
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.ChartFactory; 
import org.jfree.chart.plot.PlotOrientation; 
import org.jfree.data.xy.XYSeriesCollection; 
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.ChartUtilities;

import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;

public class FittingAffineTransforms extends ApplicationFrame {

   public FittingAffineTransforms( String applicationTitle, String chartTitle, int N ) {
      super(applicationTitle);
      JFreeChart xylineChart = ChartFactory.createXYLineChart(
         chartTitle + " Using " + N + " Points" ,
         "X Position" ,
         "Y Position" ,
         createDataset(N) ,
         PlotOrientation.VERTICAL ,
         true , true , false);
         
      ChartPanel chartPanel = new ChartPanel( xylineChart );
      chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
      final XYPlot plot = xylineChart.getXYPlot( );
      
      XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(false, true );
      renderer.setSeriesPaint( 0 , Color.BLUE );
      renderer.setSeriesPaint( 1 , Color.RED );
      renderer.setSeriesPaint( 2 , Color.GREEN );
      renderer.setSeriesPaint( 3 , Color.YELLOW );
      renderer.setSeriesStroke( 0 , new BasicStroke( 1.0f ) );
      renderer.setSeriesStroke( 1 , new BasicStroke( 1.0f ) );
      renderer.setSeriesStroke( 2 , new BasicStroke( 1.0f ) );
      renderer.setSeriesStroke( 3 , new BasicStroke( 1.0f ) );
      renderer.setSeriesShapesFilled(0, false);
      renderer.setSeriesShapesFilled(1, false);
      renderer.setSeriesShapesFilled(2, false);
      renderer.setSeriesShapesFilled(3, false);
      plot.setRenderer( renderer ); 
      setContentPane( chartPanel ); 
   }
   
   RealMatrix compose2dAffineTransform(double xScale, double yScale, double shear, double rotAngle, double xTrans, double yTrans) {
       RealMatrix rot = MatrixUtils.createRealMatrix(new double[][] {{Math.cos(rotAngle), -Math.sin(rotAngle)},
                                                                     {Math.sin(rotAngle),  Math.cos(rotAngle)}});
       
       RealMatrix sh = MatrixUtils.createRealIdentityMatrix(2);
       sh.setEntry(0, 1, shear);
       
       RealMatrix sc = MatrixUtils.createRealMatrix(2, 2);
       sc.setEntry(0, 0, xScale);
       sc.setEntry(1, 1, yScale);
       
       RealMatrix linearTransform = rot.multiply(sh).multiply(sc);
       
       RealMatrix ret = MatrixUtils.createRealMatrix(3, 3);
       ret.setSubMatrix(linearTransform.getData(), 0, 0);
       ret.setEntry(0, 2, xTrans);
       ret.setEntry(1, 2, yTrans);
       ret.setEntry(2, 0, 1.0);
       ret.setEntry(2, 1, 1.0);
       
       return ret;
   }
   
   RealMatrix toHomogeneousCoordinates(RealMatrix A) {
       int dim = A.getRowDimension();
       int N = A.getColumnDimension();
       RealMatrix ret = MatrixUtils.createRealMatrix(dim+1,N);
       RealMatrix ones = MatrixUtils.createRealMatrix(1, N).scalarAdd(1.0);
       ret.setSubMatrix(A.getData(), 0, 0);
       ret.setSubMatrix(ones.getData(), dim, 0);
       return ret;
   }
   
   RealMatrix fromHomogeneousCoordinates(RealMatrix A) {
       int dim = A.getRowDimension();
       int N = A.getColumnDimension();
       return A.getSubMatrix(0, dim-2, 0, N-1);
   }
   
   RealMatrix computeCentroid(RealMatrix A) {
       int dim = A.getRowDimension();
       int N = A.getColumnDimension();
       RealMatrix ret = MatrixUtils.createRealMatrix(dim,1);
       for (int j=0; j<N; j++) {
           ret = ret.add(A.getColumnMatrix(j));
       }
       return ret.scalarMultiply(1.0/N);
   }
   
   RealMatrix addToAllColumns(RealMatrix A, RealMatrix colMatrix ) {
       int dim = A.getRowDimension();
       int N = A.getColumnDimension();
       RealMatrix ret = MatrixUtils.createRealMatrix(dim,N);
       for (int j=0; j<N; j++) {
           ret.setColumnMatrix(j, A.getColumnMatrix(j).add(colMatrix));
       }
       return ret;
   }
   
   RealMatrix subtractFromAllColumns(RealMatrix A, RealMatrix colMatrix ) {
       int dim = A.getRowDimension();
       int N = A.getColumnDimension();
       RealMatrix ret = MatrixUtils.createRealMatrix(dim,N);
       for (int j=0; j<N; j++) {
           ret.setColumnMatrix(j, A.getColumnMatrix(j).subtract(colMatrix));
       }
       return ret;
   }
   
   RealMatrix compute3PointAffineTransform(RealMatrix A, RealMatrix B) {
       RealMatrix A_ = MatrixUtils.createRealMatrix(3, 3).scalarAdd(1.0);
       A_.setSubMatrix(A.getSubMatrix(0, 1, 0, 2).getData(), 0, 0);
       RealMatrix B_ = B.getSubMatrix(0, 1, 0, 2);
       RealMatrix inverse = new LUDecomposition(A_).getSolver().getInverse();
       RealMatrix transform = B_.multiply(inverse);
       RealMatrix ret = MatrixUtils.createRealMatrix(3, 3).scalarAdd(1.0);
       ret.setSubMatrix(transform.getData(), 0, 0);
       ret.setEntry(2, 2, 0.0);
       return ret;
   }
   
   RealMatrix computeBestAffineTransform(RealMatrix A, RealMatrix B) {
       int dim = A.getRowDimension();
       RealMatrix Acentroid = computeCentroid(A);
       RealMatrix Bcentroid = computeCentroid(B);
       RealMatrix Ahat = subtractFromAllColumns(A, Acentroid);
       RealMatrix Bhat = subtractFromAllColumns(B, Bcentroid);
       RealMatrix pinv = new SingularValueDecomposition(Ahat.multiply(Ahat.transpose())).getSolver().getInverse();
       RealMatrix linearTransform = Bhat.multiply(Ahat.transpose()).multiply(pinv);
       RealMatrix translation = Bcentroid.subtract(linearTransform.multiply(Acentroid));
       RealMatrix ret = MatrixUtils.createRealMatrix(dim+1, dim+1).scalarAdd(1.0);
       ret.setSubMatrix(linearTransform.getData(), 0, 0);
       ret.setSubMatrix(translation.getData(), 0, dim);
       ret.setEntry(dim, dim, 0.0);
       return ret;
   }

   RealMatrix applyAffineTransform(RealMatrix affineTransform, RealMatrix data) {
       return fromHomogeneousCoordinates(affineTransform.multiply(toHomogeneousCoordinates(data)));
   }
   
   private XYDataset createDataset(int N) {
      double measurementError = 5.0; //size of measurement errors
      
      //Create a random set of points
      double[][] src = new double[2][N];
      for (int i=0; i<2; i++) {
          for (int j=0; j<N; j++) {
              src[i][j] = 100.0*Math.random()-50.0;
          }
      }
      RealMatrix source = MatrixUtils.createRealMatrix(src);
      
      //Create a random Affine transform
      double xScale = 1.0 + 0.2 * (Math.random()-0.5);
      double yScale = 1.0 + 0.2 * (Math.random()-0.5);
      double shear = 0.0 + 0.2 * (Math.random()-0.5);
      double rotAngle = 2.0 * Math.PI * Math.random();
      double xTrans = 200 * (Math.random()-0.5);
      double yTrans = 200 * (Math.random()-0.5);
      RealMatrix trueAT = compose2dAffineTransform(xScale, yScale, shear, rotAngle, xTrans, yTrans);
      
      RealMatrix trueData = applyAffineTransform(trueAT, source);
      
      for (int i=0; i<2; i++) {
          for (int j=0; j<N; j++) {
              src[i][j] = measurementError*2.0*(Math.random()-0.5);
          }
      }
      RealMatrix measurementErr = MatrixUtils.createRealMatrix(src);
      
      RealMatrix measuredData = trueData.add(measurementErr);
      
      RealMatrix recoveredAT = computeBestAffineTransform(source, measuredData);
      RealMatrix originalAT = compute3PointAffineTransform(source, measuredData);
      
      RealMatrix recoveredData = applyAffineTransform(recoveredAT, source);
      RealMatrix originalData = applyAffineTransform(originalAT, source);

      RealMatrix recoveredErr = recoveredData.subtract(trueData);
      double recoveredRMS = Math.sqrt(recoveredErr.transpose().multiply(recoveredErr).getTrace() / N);
      
      RealMatrix originalErr = originalData.subtract(trueData);
      double originalRMS = Math.sqrt(originalErr.transpose().multiply(originalErr).getTrace() / N);
      
      System.out.println("New RMS error = " + recoveredRMS + ", original = " + originalRMS);
      
      final XYSeries truth = new XYSeries( "Truth" );
      final XYSeries measured = new XYSeries( "Measured" );
      final XYSeries recovered = new XYSeries( "Recovered" );
      final XYSeries original = new XYSeries( "Original" );
      for (int j=0; j<N; j++) {
          truth.add(trueData.getEntry(0, j), trueData.getEntry(1, j));
          measured.add(measuredData.getEntry(0, j), measuredData.getEntry(1, j));
          recovered.add(recoveredData.getEntry(0, j), recoveredData.getEntry(1, j));
          original.add(originalData.getEntry(0, j), originalData.getEntry(1, j));
      }
      
      final XYSeriesCollection dataset = new XYSeriesCollection( );          
      dataset.addSeries( truth );          
      dataset.addSeries( measured );          
      dataset.addSeries( recovered );
      dataset.addSeries( original );
      return dataset;
   }

   public static void main( String[ ] args ) {
       FittingAffineTransforms chart = new FittingAffineTransforms("Demo of Fitting Affine Transforms",
         "Best Affine Transform", 4);
      chart.pack( );          
      RefineryUtilities.centerFrameOnScreen( chart );          
      chart.setVisible( true ); 
   }
}