public class Edge {
	private Node fromNode;
	private Node toNode;

	public Node getFromNode() {
		return fromNode;
	}

	public Node getToNode() {
		return toNode;
	}

	public void setFromNode(Node fromNode) {
		this.fromNode = fromNode;
	}

	public void setToNode(Node toNode) {
		this.toNode = toNode;
	}


	@Override
	public String toString() {
		return String.format("%s->%s", fromNode, toNode);
	}

	@Override
	public int hashCode() {
		return (fromNode.toString() + toNode.toString()).hashCode();
	}

	@Override
	public boolean equals(Object other) {
		if (other instanceof Edge) {
			Edge o = (Edge) other;
			return fromNode.toString().equals(o.fromNode.toString()) && toNode.toString().equals(o.toNode.toString());
		}
		else
			return false;
	}
}
