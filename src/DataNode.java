import java.util.ArrayList;
import java.util.List;


public class DataNode implements Node {
	private String data;
	private List<Edge> incomingEdgeList = new ArrayList<Edge>();
	private List<Edge> outgoingEdgeList = new ArrayList<Edge>();
	private int layer = -1;

	public DataNode(String data) {
		this.data = data;
	}

	public List<Edge> getIncomingEdgeList() {
		return incomingEdgeList;
	}

	public List<Edge> getOutgoingEdgeList() {
		return outgoingEdgeList;
	}

	public String getData() {
		return data;
	}

	public int getLayer() {
		return layer;
	}

	public void setLayer(int layer) {
		if (this.layer == -1 || layer < this.layer)
			this.layer = layer;
	}

	@Override
	public String toString() {
		return data;
	}

	@Override
	public int hashCode() {
		return data.hashCode();
	}

	@Override
	public boolean equals(Object other) {
		if (other instanceof DataNode) {
			DataNode o = (DataNode) other;
			return data.equals(o.data);
		}
		else
			return false;
	}
}
