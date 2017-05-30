package org.gephi.spectralspatialisation;

import org.gephi.graph.api.Column;
import org.gephi.graph.api.Edge;
import org.gephi.graph.api.GraphModel;
import org.gephi.graph.api.Graph;
import org.gephi.graph.api.Node;

import org.gephi.layout.spi.Layout;
import org.gephi.layout.spi.LayoutBuilder;
import org.gephi.layout.spi.LayoutProperty;

import org.openide.util.NbBundle;
import org.openide.util.Exceptions;

import java.util.ArrayList;
import java.util.List;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;

/**
* Implementing a method detailed by Jerome Kunegis in the linked PDF files
* Spatialisation relying on a Graph Laplacian built on adjacency matrix
* @author CERI
* @see http://epubs.siam.org/doi/pdf/10.1137/1.9781611972801.49
* @see https://arxiv.org/pdf/1402.6865.pdf
*/
public class SpectralSpatialisation implements Layout {
	private final LayoutBuilder build;
	protected GraphModel model;
	protected Graph graph;

	private final double areaFX  = 512.0;

	private double _scale            = 1.0;
	private double _thirdComponent   = 0.0;
	private double _fourthComponent  = 0.0;

	private DoubleMatrix2D ideal     =     null;

	private String _weight_attribute = "weight";
	private boolean _weight_ignore   =    false;
	private String _sign_attribute   =   "sign";
	private boolean _sign_ignore     =    false;
	public static String[] rows      = {"weight","sign"};

	private boolean ran = false;

	public SpectralSpatialisation(LayoutBuilder build) {
		this.build = build;
	}

	/**
	 * This method sets the variables back to their initial values.
	 */
	@Override
	public void resetPropertiesValues() {
		_scale = 1.0;
		_thirdComponent = 0.0;

		// Reset rows
		rows = new String[model.getEdgeTable().countColumns()];
		int i = 0;
		// Fill rows with correct names
		for (Column c : model.getEdgeTable()) {rows[i]=c.getTitle(); i++;}
	}
	/**
	 * This method initializes the variables that will be used to place the nodes. (aka process the data)
	 */
	@Override
	public void initAlgo() {
		this.graph = model.getGraphVisible();
		graph.readLock();
		
		try {
			Node[] nodes = graph.getNodes().toArray();
			Edge[] edges = graph.getEdges().toArray();

			// PROCESSING

			// Create the Adjacency matrix (The sign*weight of the {u~v} edge or zero if no edge)
			// Create the Diagonal matrix (Number of neighbors) (in fact not)
			// L = D-A
			// That's the theory.

			// Practical :
			// Generate -A
			// Manually compute D
			// Manually add D to -A, resulting in L.

			// Get eigenvalues
			// Get eigenvectors matching the n least non-zero eigenvalues
			// take them as X = xy(zw) (Ideal position)

			// Four components are taken, the third because of known problems in case of unsigned graph
			// And the fourth for fun, doesn't cost much more.

			DoubleMatrix2D laplacian;
			laplacian = new SparseDoubleMatrix2D(nodes.length,nodes.length);
			if (nodes.length * nodes.length > 4 * edges.length) {
				laplacian = new DenseDoubleMatrix2D(nodes.length,nodes.length);
			}
			double[] diag = new double[nodes.length];
			for (Edge e : edges) {
				double sign = getSign(e);
				double weight = getWeight(e);
				double val = sign*weight;
				laplacian.set(e.getSource().getStoreId(), e.getTarget().getStoreId(), -val);
				laplacian.set(e.getTarget().getStoreId(), e.getSource().getStoreId(), -val);
				if (val < 0) {val = -val;}
				// Manual computation of D, Sum[absoulte values of neighboring edges]
				diag[e.getSource().getStoreId()] += val;
				diag[e.getTarget().getStoreId()] += val;
			}
			for (Node n : nodes) {
				int id = n.getStoreId();
				double w = laplacian.get(id, id);
				laplacian.set(id, id, w + diag[id]);
			}

			// Do the decomposition, giving both values and vectors
			EigenvalueDecomposition eigen = new EigenvalueDecomposition(laplacian);
			// First get the values in order to keep only the n least (by absolute value, and non-zero)
			DoubleMatrix2D values = eigen.getD();

			int limit = 4; // How much do we keep
			int[] least_index = new int[limit];       // Sorting memory stuff
			double[] least_value = new double[limit]; // Sorting memory stuff
			for (int i = 0 ; i < limit ; i++) {least_index[i] = -1;}

			int index;
			double value;

			for (int i = 0 ; i < nodes.length ; i ++) {
				index = i;
				value = values.get(i,i);
				// We need absolute value
				if (value < 0) {value = -value;}
				// If it's the null value, skip it
				if (value == 0) {}
				// Else, it becomes more interesting
				else {
					for(int k = 0 ; k < limit ; k++) {
						// If this index is not yet set
						if (least_index[k] == -1) {
							least_index[k] = index;
							least_value[k] = value;
							// Here it is important to break
							// As indexes are always sorted by best-values (at any step)
							// then the next index is not yet set.
							// So the values that would be given to the next loop (same as the ones of this iteration)
							// would be accepted by the next index, leading to an eigen value being selected multiple times (not acceptable)
							// It would not be a problem if the indexes were changed enough across successive loops, they would end up
							// holding different eigen values.
							// But let's say the very best non-zero values comes first, then all indexes would be filled will it by the end of the
							// first loop and non of the other eigen values would fit in any least_index after that
							break;
						}
						// If it's set but as a less interesting eigen valeus associated with it
						else if (least_value[k] > value) {
							// Keep current values, as they will be passed to next iteration
							int si = least_index[k];
							double sv = least_value[k];
							// Set new best values
							least_index[k] = index;
							least_value[k] = value;
							// Here it is important to keep on looping
							// If an index is changed, let's say least_index[0]
							// Then as indexes were already sorted by best-values
							// least_index[1] must store an index worse than current least_index[0]
							// So if an index is changed, it has to push its value into the next one
							// Cascading until first unset least_index
							index = si; value = sv;
						}
					}
				}
			}
			// Now we got to select vectors
			DoubleMatrix2D vectors = eigen.getV();
			ideal = vectors.viewSelection(null,least_index);
		}
		finally {
			graph.readUnlockAll();
		}
	}
	/**
	 * Based on initAlgo(), this method makes the nodes reach their ideal position
	 * Third component can be used to alter to behavior.
	 */
	@Override
	public void goAlgo() {
		this.graph = model.getGraphVisible();
		graph.readLock();
		try {
			float fx = (float)(areaFX * _scale);
			Node[] nodes = graph.getNodes().toArray();
			for (Node n : nodes) {
				// Place each node where it belongs
				n.setX( fx * (float)(ideal.get(n.getStoreId(), 0) + _thirdComponent  * ideal.get(n.getStoreId(), 2)));
				n.setY(-fx * (float)(ideal.get(n.getStoreId(), 1) + _fourthComponent * ideal.get(n.getStoreId(), 3)));
			}
		}
		finally {
			graph.readUnlockAll();
			ran = true;
		}
	}

	@Override
	public void endAlgo() {
		graph.readLock();
		try {
			for (Node n : graph.getNodes()) {
				n.setLayoutData(null);
			}
		}
		finally {
			graph.readUnlockAll();
			ran = false;
		}
	}

	@Override
	public boolean canAlgo() {
		return !ran;
	}

	@Override
	public LayoutProperty[] getProperties() {
		List<LayoutProperty> properties = new ArrayList<>();
		final String SPECTRALSPATIALISATION = "Spectral Spatialisation (Signed Graph Laplacian)";

		try {
			properties.add(LayoutProperty.createProperty(
				this, String.class,
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.weight_attribute.name"),
				SPECTRALSPATIALISATION,
				"SpectralSpatialisation.weight_attribute.name",
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.weight_attribute.desc"),
				"getWeightAttribute", "setWeightAttribute", CustomComboBoxEditor.class));
			properties.add(LayoutProperty.createProperty(
				this, boolean.class,
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.weight_ignore.name"),
				SPECTRALSPATIALISATION,
				"SpectralSpatialisation.weight_ignore.name",
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.weight_ignore.desc"),
				"getWeightIgnore", "setWeightIgnore"));

			properties.add(LayoutProperty.createProperty(
				this, String.class,
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.sign_attribute.name"),
				SPECTRALSPATIALISATION,
				"SpectralSpatialisation.sign_attribute.name",
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.sign_attribute.desc"),
				"getSignAttribute", "setSignAttribute", CustomComboBoxEditor.class));
			properties.add(LayoutProperty.createProperty(
				this, boolean.class,
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.sign_ignore.name"),
				SPECTRALSPATIALISATION,
				"SpectralSpatialisation.sign_ignore.name",
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.sign_ignore.desc"),
				"getSignIgnore", "setSignIgnore"));

			properties.add(LayoutProperty.createProperty(
				this, double.class,
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.scale.name"),
				SPECTRALSPATIALISATION,
				"SpectralSpatialisation.scale.name",
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.scale.desc"),
				"getScale", "setScale"));

			properties.add(LayoutProperty.createProperty(
				this, double.class,
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.thirdComponent.name"),
				SPECTRALSPATIALISATION,
				"SpectralSpatialisation.thirdComponent.name",
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.thirdComponent.desc"),
				"getThirdComponent", "setThirdComponent"));
			properties.add(LayoutProperty.createProperty(
				this, double.class,
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.fourthComponent.name"),
				SPECTRALSPATIALISATION,
				"SpectralSpatialisation.fourthComponent.name",
				NbBundle.getMessage(SpectralSpatialisation.class, "SpectralSpatialisation.fourthComponent.desc"),
				"getFourthComponent", "setFourthComponent"));
			
		} catch (Exception e) {
			Exceptions.printStackTrace(e);
		}
		return properties.toArray(new LayoutProperty[0]);
	}

	public String getWeightAttribute() {
		return _weight_attribute;
	}
	public void setWeightAttribute(String attribute) {
		_weight_attribute = attribute;
	}
	private double getWeight(Edge e) {
		if (_weight_ignore) {return 1.0;}
		Object attr = e.getAttribute(_weight_attribute);
		if (attr instanceof Boolean) {
			return (((boolean)attr) ? 1.0 : 0.0);
		}
		if (attr instanceof Integer) {return Math.abs((double)((int)attr));}
		if (attr instanceof Float)   {return Math.abs((double)((float)attr));}
		if (attr instanceof Double)  {return Math.abs((double)attr);}
		return 1.0;
	}
	public boolean getWeightIgnore() {
		return _weight_ignore;
	}
	public void setWeightIgnore(boolean ignore) {
		_weight_ignore = ignore;
	}

	public String getSignAttribute() {
		return _sign_attribute;
	}
	public void setSignAttribute(String attribute) {
		_sign_attribute = attribute;
	}
	private double getSign(Edge e) {
		if (_sign_ignore) {return 1.0;}
		Object attr = e.getAttribute(_sign_attribute);
		if (attr instanceof Boolean) {
			return (double)(((boolean)attr) ? 1.0 : -1.0);
		}
		if (attr instanceof Integer) {return ((((int)attr)    >= 0) ? 1.0 : -1.0);}
		if (attr instanceof Float)   {return ((((float)attr)  >= 0) ? 1.0 : -1.0);}
		if (attr instanceof Double)  {return ((((double)attr) >= 0) ? 1.0 : -1.0);}
		return 1.0;
	}
	public boolean getSignIgnore() {
		return _sign_ignore;
	}
	public void setSignIgnore(boolean ignore) {
		_sign_ignore = ignore;
	}

	public double getScale() {
		return _scale;
	}
	public void setScale(double scale) {
		if (scale <= 0.0) {_scale = 1.0;}
		else {_scale = scale;}
	}
	public double getThirdComponent() {
		return _thirdComponent;
	}
	public void setThirdComponent(double thirdComponent) {
		_thirdComponent = thirdComponent;
	}

	public double getFourthComponent() {
		return _fourthComponent;
	}
	public void setFourthComponent(double fourthComponent) {
		_fourthComponent = fourthComponent;
	}


	@Override
	public LayoutBuilder getBuilder() {
		return build;
	}
	@Override
	public void setGraphModel(GraphModel model) {
		this.model = model;
	}
}