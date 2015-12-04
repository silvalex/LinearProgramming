import gurobi.*;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.w3c.dom.Text;
import org.xml.sax.SAXException;


/**
 * This Web service composition strategy was implemented as a recreation of the following paper:
 *
 * Gabrel, Virginie, Maude Manouvrier, and Cécile Murat. "Optimal and Automatic Transactional Web
 * Service Composition with Dependency Graph and 0-1 Linear Programming." Service-Oriented
 * Computing. Springer Berlin Heidelberg, 2014. 108-122.
 *
 * @author sawczualex
 */
public class WebServiceComposition {
	// Constants with of order of QoS attributes
	public static final int TIME = 0;
	public static final int COST = 1;
	public static final int AVAILABILITY = 2;
	public static final int RELIABILITY = 3;

	String dotFileName = "wsc";
	String imageFileName = "wsc";

	public final double minAvailability = 0.0;
	public double maxAvailability = -1.0;
	public final double minReliability = 0.0;
	public double maxReliability = -1.0;
	public double minTime = Double.MAX_VALUE;
	public double maxTime = -1.0;
	public double minCost = Double.MAX_VALUE;
	public double maxCost = -1.0;

	public double w1 = 0.25;
	public double w2 = 0.25;
	public double w3 = 0.25;
	public double w4 = 0.25;

	public Map<String, ServiceNode> serviceMap = new HashMap<String, ServiceNode>();
	public Set<ServiceNode> relevantServices = new HashSet<ServiceNode>();
	public Map<String, DataNode> relevantData = new HashMap<String, DataNode>();
	public Map<String, TaxonomyNode> taxonomyMap = new HashMap<String, TaxonomyNode>();
	public Set<String> taskInput;
	public Set<String> dependencyGraphInputs = new HashSet<String>();
	public Set<String> taskOutput;

	public Map<String, ServiceNode> dependencyGraphNodes = new HashMap<String, ServiceNode>();
	public Map<String, DataNode> dependencyGraphData = new HashMap<String, DataNode>();
	public Map<String, Edge> dependencyGraphEdges = new HashMap<String, Edge>();
	public Map<String, GRBVar> variables = new HashMap<String, GRBVar>();

	public ServiceNode startNode;


	public static void main(String[] args) {
	    if (args.length < 3) {
	        System.out.println("Usage: java WebServiceComposition problem.xml services-output.xml taxonomy.xml");
	        System.exit(1);
	    }
	    new WebServiceComposition(args[0], args[1], args[2]);
	}

	public WebServiceComposition(String problem, String servicesOutput, String taxonomy) {
		// Parse files
		parseWSCTaskFile(problem);
		parseWSCServiceFile(servicesOutput);
		parseWSCTaxonomyFile(taxonomy);

		// Find concepts for instances
		findConceptsForInstances();

		// Populate the taxonomy tree with services
		populateTaxonomyTree();

		// Find relevant nodes for the dependency graph (both services and data)
		getRelevantNodes(serviceMap, relevantServices, relevantData, taskInput, taskOutput);

		// Calculate boundaries for normalisation
		calculateNormalisationBounds(relevantServices);

		// Create dependency graph
		System.out.println("Dependency graph inputs:");
		System.out.println(dependencyGraphInputs);
		System.out.println("Dependency graph:");
		createDependencyGraph(relevantServices, relevantData, dependencyGraphEdges, dependencyGraphNodes, dependencyGraphData);
		System.out.println(edgesToDot(dependencyGraphEdges));

		try {
			GRBEnv env = new GRBEnv();
			GRBModel model = new GRBModel(env);

			// Create model variables
			for (ServiceNode s : dependencyGraphNodes.values())
				variables.put(s.toString(), model.addVar(0.0, 1.0, 0.0, GRB.BINARY, s.toString()));
			for (Edge e : dependencyGraphEdges.values())
				variables.put(e.toString(), model.addVar(0.0, 1.0, 0.0, GRB.BINARY, e.toString()));

			// Integrate variables into model
		      model.update();

		      // Add constraints
		      GRBLinExpr leftSide;
		      GRBVar[] vars;
		      GRBLinExpr rightSide;

		      /* (C1) For an output edge of a service to be used (i.e. included in the composition),
		       * all inputs edges of that service must also be included in the composition (i.e. all service
		       * inputs must be fulfilled for its output to be usable)
		       */
		      for (ServiceNode sn : dependencyGraphNodes.values()) {
    			  // Predecessor edges
    			  leftSide = new GRBLinExpr();
    			  List<Edge> predecessorEdges = sn.getIncomingEdgeList();
    			  double numPredecessors = predecessorEdges.size();
    			  vars = new GRBVar[(int)numPredecessors];
    			  for (int i = 0; i < numPredecessors; i++)
    				  vars[i] = variables.get(predecessorEdges.get(i).toString());
    			  leftSide.addTerms(null, vars);

	    		  // For all successor edges
	    		  for (Edge e : sn.getOutgoingEdgeList()) {
	    			  rightSide = new GRBLinExpr();
	    			  rightSide.addTerm(numPredecessors, variables.get(e.toString()));
	    			  String constName = "C1_" + sn.toString() + "_" + e.toString();
			    	  model.addConstr(leftSide, GRB.GREATER_EQUAL, rightSide, constName);
	    		  }
		      }

		      /* (C2) For the outgoing edges of a data node to be used, at least one incoming edge to that
		       * data node must also be included in the composition (excluding data provided by user).
		       */
		      // For all data nodes not given as input by the user
		      for (DataNode dn : dependencyGraphData.values()) {
		    	  if (!dependencyGraphInputs.contains(dn.toString())) {

	    			  // Predecessor edges
	    			  leftSide = new GRBLinExpr();
	    			  List<Edge> predecessorEdges = dn.getIncomingEdgeList();
	    			  vars = new GRBVar[predecessorEdges.size()];
	    			  for (int i = 0; i < predecessorEdges.size(); i++)
	    				  vars[i] = variables.get(predecessorEdges.get(i).toString());
	    			  leftSide.addTerms(null, vars);

		    		  // For all successor edges
		    		  for (Edge e : dn.getOutgoingEdgeList()) {
		    			  rightSide = new GRBLinExpr();
		    			  rightSide.addTerm(1.0, variables.get(e.toString()));
		    			  String constName = "C2_" + dn.toString() + "_" + e.toString();
				    	  model.addConstr(leftSide, GRB.GREATER_EQUAL, rightSide, constName);
		    		  }
		    	  }
		      }

		      /* (C3) For all user output data nodes, at least one incoming edge should be included
		       * in the composition.
		       */
		      for (DataNode dn : dependencyGraphData.values()) {
		    	  if (taskOutput.contains(dn.toString())) {
		    		  leftSide = new GRBLinExpr();
		    		  List<Edge> predecessorEdges = dn.getIncomingEdgeList();
		    		  vars = new GRBVar[predecessorEdges.size()];
	    			  for (int i = 0; i < predecessorEdges.size(); i++)
	    				  vars[i] = variables.get(predecessorEdges.get(i).toString());
	    			  leftSide.addTerms(null, vars);
	    			  String constName = "C3_" + dn.toString();
			    	  model.addConstr(leftSide, GRB.GREATER_EQUAL, 1.0, constName);
		    	  }
		      }

		      /* (C4) The outgoing edges of a service can only be included in a composition if the service
		       *  itself is also included.
		       */
		      for (ServiceNode sn : dependencyGraphNodes.values()) {
		    	  leftSide = new GRBLinExpr();
		    	  List<Edge> successorEdges = sn.getOutgoingEdgeList();
		    	  double numSuccessors = successorEdges.size();
		    	  vars = new GRBVar[(int)numSuccessors];
    			  for (int i = 0; i < numSuccessors; i++)
    				  vars[i] = variables.get(successorEdges.get(i).toString());
    			  leftSide.addTerms(null, vars);

    			  rightSide = new GRBLinExpr();
    			  rightSide.addTerm(numSuccessors, variables.get(sn.toString()));
    			  String constName = "C4_" + sn.toString();
		    	  model.addConstr(leftSide, GRB.LESS_EQUAL, rightSide, constName);
		      }

		      /* (C5) For all edges included in the composition, the edge must flow from a node with a
		       * smaller topology score to a node with a larger topology score (cycle prevention).
		       */
		      for (Edge e : dependencyGraphEdges.values()) {
		    	  double topologyDiff = e.getToNode().getLayer() - e.getFromNode().getLayer();
		    	  double vertexSize = dependencyGraphNodes.size() + dependencyGraphData.size();

		    	  rightSide = new GRBLinExpr();// 1 -|X| + |X| * x_{i,j}
		    	  rightSide.addConstant(1.0);
		    	  rightSide.addConstant(-vertexSize);
		    	  rightSide.addTerm(vertexSize, variables.get(e.toString()));

		    	  String constName = "C5_" + e.toString();
		    	  model.addConstr(topologyDiff, GRB.GREATER_EQUAL, rightSide, constName);
		      }

		      /* (C6) The topology score of all user inputs is 0.
		       */
		      for (DataNode dn : dependencyGraphData.values()) {
		    	  if (dependencyGraphInputs.contains(dn.toString())) {
		    		  leftSide = new GRBLinExpr();
		    		  leftSide.addConstant(dn.getLayer());
		    		  String constName = "C6_" + dn.toString();
			    	  model.addConstr(leftSide, GRB.EQUAL, 0, constName);
		    	  }
		      }

		      // Update model
		      model.update();

		      /* (OBJECTIVE) Minimize according to the sum of scores from each individual
		       * service included in the composition. The score is calculated as a weighted sum of
		       * its normalized quality aspects, without taking any workflow constructs into consideration.
		       */
		      GRBLinExpr obj = new GRBLinExpr();

		      for (ServiceNode sn : dependencyGraphNodes.values()) {
		    	  obj.addTerm(calculateQos(sn), variables.get(sn.toString()));
		      }
		      model.setObjective(obj, GRB.MINIMIZE);

		      // Update model
		      model.update();

		      // Optimize model
		      model.optimize();

		      // Write model to file
		      model.write("wsc.lp");

		      // Generate graph solution using dot syntax
		      Map<String, Double> valueMap = new HashMap<String, Double>();

		      GRBVar[] varArray = new GRBVar[variables.size()];
		      String[] varNames = new String[variables.size()];
		      variables.values().toArray(varArray);
		      variables.keySet().toArray(varNames);
		      double[] values = model.get(GRB.DoubleAttr.X, varArray);

		      for (int i = 0; i < values.length; i++)
		    	  valueMap.put(varNames[i], values[i]);

		      StringBuilder builder = new StringBuilder();
		      builder.append("digraph g{");

			for (Entry<String, Double> e : valueMap.entrySet()) {
				if (e.getKey().contains("->") && e.getValue() == 1.0) {
					// Add edge to the solution
					builder.append(e.getKey());
					builder.append(";\n");

					// Check that service ends of edges are also included in the
					// solution
					Edge edge = dependencyGraphEdges.get(e.getKey());
					ServiceNode serv = null;
					if (edge.getFromNode() instanceof ServiceNode) {
						serv = (ServiceNode) edge.getFromNode();
					} else if (edge.getToNode() instanceof ServiceNode) {
						serv = (ServiceNode) edge.getToNode();
					} else {
						throw new RuntimeException(String.format(
								"Edge '%s' does not link to a service node!",
								edge.toString()));
					}
					if (valueMap.get(serv.toString()) == 0.0) {
						builder.append(serv.toString());
						builder.append(" [color=\"#FF0000\"];\n");
					}
				}
			}

			builder.append("}");
			try {
				dotFileName = String.format("%s.dot", dotFileName);
				imageFileName = String.format("%s.png", imageFileName);

				FileWriter writer = new FileWriter(new File(dotFileName));
				writer.append(builder.toString());
				writer.close();
				java.lang.Runtime.getRuntime().exec(String.format("dot -Tpng %s -o %s", dotFileName, imageFileName));
				System.out.println(String.format("Solution image '%s' created.", imageFileName));
				java.lang.Runtime.getRuntime().exec(String.format("gwenview %s", imageFileName));

			} catch (IOException e1) {
				throw new RuntimeException("Problem creating solution graph!");
			}
			// Dispose of model and environment
			model.dispose();
			env.dispose();

		}
		catch (GRBException e) {
			System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
	}

	private double calculateQos(ServiceNode s) {
		double[] qos = s.getQos();
        double a = qos[AVAILABILITY];
        double r = qos[RELIABILITY];
        double t = qos[TIME];
        double c = qos[COST];

        a = normaliseAvailability(a);
        r = normaliseReliability(r);
        t = normaliseTime(t);
        c = normaliseCost(c);

		return w1 * a + w2 * r + w3 * t + w4 * c;
	}

	private double normaliseAvailability(double availability) {
		if (maxAvailability - minAvailability == 0.0)
			return 1.0;
		else
			return (maxAvailability - availability)/(maxAvailability - minAvailability);
	}

	private double normaliseReliability(double reliability) {
		if (maxReliability - minReliability == 0.0)
			return 1.0;
		else
			return (maxReliability - reliability)/(maxReliability - minReliability);
	}

	private double normaliseTime(double time) {
		if (maxTime - minTime == 0.0)
			return 1.0;
		else
			return (time - minTime)/(maxTime - minTime);
	}

	private double normaliseCost(double cost) {
		if (maxCost - minCost == 0.0)
			return 1.0;
		else
			return (cost - minCost)/(maxCost - minCost);
	}

	private void calculateNormalisationBounds(Set<ServiceNode> services) {
		for(ServiceNode service: services) {
			double[] qos = service.getQos();

			// Availability
			double availability = qos[AVAILABILITY];
			if (availability > maxAvailability)
				maxAvailability = availability;

			// Reliability
			double reliability = qos[RELIABILITY];
			if (reliability > maxReliability)
				maxReliability = reliability;

			// Time
			double time = qos[TIME];
			if (time > maxTime)
				maxTime = time;
			if (time < minTime)
				minTime = time;

			// Cost
			double cost = qos[COST];
			if (cost > maxCost)
				maxCost = cost;
			if (cost < minCost)
				minCost = cost;
		}
		// Adjust max. cost and max. time based on the number of services in shrunk repository
		maxCost *= services.size();
		maxTime *= services.size();

	}

	/**
	 * Converts input, output, and service instance values to their corresponding
	 * ontological parent.
	 */
	private void findConceptsForInstances() {
		Set<String> temp = new HashSet<String>();

		for (String s : taskInput)
			temp.add(taxonomyMap.get(s).parents.get(0).value);
		taskInput.clear();
		taskInput.addAll(temp);

		temp.clear();
		for (String s : taskOutput)
				temp.add(taxonomyMap.get(s).parents.get(0).value);
		taskOutput.clear();
		taskOutput.addAll(temp);

		for (ServiceNode s : serviceMap.values()) {
			temp.clear();
			Set<String> inputs = s.getInputs();
			for (String i : inputs)
				temp.add(taxonomyMap.get(i).parents.get(0).value);
			inputs.clear();
			inputs.addAll(temp);

			temp.clear();
			Set<String> outputs = s.getOutputs();
			for (String o : outputs)
				temp.add(taxonomyMap.get(o).parents.get(0).value);
			outputs.clear();
			outputs.addAll(temp);
		}
	}

	/**
	 * For the services provided, creates a graph that connects all the possible
	 * input and output pairs between any two nodes.
	 *
	 * @param services
	 * @return graph of input-output matches
	 */
	private void createDependencyGraph(Set<ServiceNode> relevantServices, Map<String, DataNode> relevantData, Map<String, Edge> edgeMap, Map<String, ServiceNode> serviceMap, Map<String, DataNode> dataMap) {
		// Populate dependency graph inputs by querying the taxonomy
		for (String i : taskInput) {
			for (Set<String> s : taxonomyMap.get(i).servicesWithInput.values()) {
				dependencyGraphInputs.addAll(s);
			}
		}
		
		// For all services, add their inputs as data nodes (with edges to services)
		for (ServiceNode serv : relevantServices) {
			serviceMap.put(serv.toString(), serv);
			for (String i : serv.getInputs()) {
				DataNode dn = relevantData.get(i);
				dataMap.put(i, dn);
				Edge e = new Edge();
				e.setFromNode(dn);
				e.setToNode(serv);
				dn.getOutgoingEdgeList().add(e);
				serv.getIncomingEdgeList().add(e);
				edgeMap.put(e.toString(), e);
			}
		}
		// Also add output data nodes without making connections
		for (String o : taskOutput) {
			DataNode dn = relevantData.get(o);
			dataMap.put(o, dn);
		}

		// For all data nodes (except start ones), find services that satisfy them and create edges
		for (DataNode dn : dataMap.values()) {
			if (!dependencyGraphInputs.contains(dn.toString())) {
				List<ServiceNode> services = taxonomyMap.get(dn.toString()).servicesWithOutput;
				for (ServiceNode serv : services) {
					if (relevantServices.contains(serv)) {
						Edge e = new Edge();
						e.setFromNode(serv);
						e.setToNode(dn);
						serv.getOutgoingEdgeList().add(e);
						dn.getIncomingEdgeList().add(e);
						edgeMap.put(e.toString(), e);
					}
				}
			}
		}
		// Ensure that all dependency graph inputs have layer with topology of 0
		for (String i : dependencyGraphInputs) {
			dataMap.get(i).setLayer(0);
		}
	}

	/**
	 * Goes through the service list and retrieves only those services which
	 * could be part of the composition task requested by the user.
	 *
	 * @param serviceMap
	 * @return relevant services
	 */
	private void getRelevantNodes(Map<String,ServiceNode> serviceMap, Set<ServiceNode> relevantServices, Map<String,DataNode> relevantData, Set<String> inputs, Set<String> outputs) {
		// Copy service map values to retain original
		Set<ServiceNode> services = new HashSet<ServiceNode>(serviceMap.values());

		int layer = 1;
		Set<String> cSearch = new HashSet<String>(inputs);
		Set<ServiceNode> sSet = new HashSet<ServiceNode>();
		Set<ServiceNode> sFound = discoverService(services, cSearch, layer);
		boolean firstServiceLayer = true;

		while (!sFound.isEmpty()) {
			sSet.addAll(sFound);
			services.removeAll(sFound);
			for (ServiceNode s: sFound) {
				cSearch.addAll(s.getOutputs());
				if (firstServiceLayer) {
					//dependencyGraphInputs.addAll(s.getInputs()); XXX

					// Add start node to the taxonomy
					//startNode = new ServiceNode("start", null, null, dependencyGraphInputs); XXX
					//populateOutputs(startNode); XXX
				}
				for (String i : s.getInputs()) {
					DataNode dn = new DataNode(i);
					dn.setLayer(layer-1);
					relevantData.put(i,dn);
				}
			}
			layer+=2;
			sFound.clear();
			firstServiceLayer = false;
			sFound = discoverService(services, cSearch, layer);
		}

		if (isSubsumed(outputs, cSearch)) {
			relevantServices.addAll(sSet);
			for (String o : outputs) {
				DataNode dn = new DataNode(o);
				dn.setLayer(layer-1);
				relevantData.put(o,dn);
			}
		}
		else {
			throw new RuntimeException("It is impossible to perform a composition using the services and settings provided.");
		}
	}

	/**
	 * Discovers all services from the provided collection whose
	 * input can be satisfied either (a) by the input provided in
	 * searchSet or (b) by the output of services whose input is
	 * satisfied by searchSet (or a combination of (a) and (b)).
	 *
	 * @param services
	 * @param searchSet
	 * @return set of discovered services
	 */
	private Set<ServiceNode> discoverService(Collection<ServiceNode> services, Set<String> searchSet, int layer) {
		Set<ServiceNode> found = new HashSet<ServiceNode>();
		for (ServiceNode s: services) {
			if (isSubsumed(s.getInputs(), searchSet))
				found.add(s);
				s.setLayer(layer);
		}
		return found;
	}

	/**
	 * Checks whether set of inputs can be completely satisfied by the search
	 * set, making sure to check descendants of input concepts for the subsumption.
	 *
	 * @param inputs
	 * @param searchSet
	 * @return true if search set subsumed by input set, false otherwise.
	 */
	public boolean isSubsumed(Set<String> inputs, Set<String> searchSet) {
		boolean satisfied = true;
		for (String input : inputs) {
			Set<String> subsumed = taxonomyMap.get(input).getSubsumedConcepts();
			if (!isIntersection( searchSet, subsumed )) {
				satisfied = false;
				break;
			}
		}
		return satisfied;
	}

    private static boolean isIntersection( Set<String> a, Set<String> b ) {
        for ( String v1 : a ) {
            if ( b.contains( v1 ) ) {
                return true;
            }
        }
        return false;
    }

	/**
	 * Populates the taxonomy tree by associating services to the
	 * nodes in the tree.
	 */
	private void populateTaxonomyTree() {
		for (ServiceNode s: serviceMap.values()) {
			addServiceToTaxonomyTree(s);
		}
	}

	private void addServiceToTaxonomyTree(ServiceNode s) {
		// Populate outputs
		populateOutputs(s);
		// Populate inputs
		populateInputs(s);
	}

	private void populateOutputs(ServiceNode s) {
		// Populate outputs
	    Set<TaxonomyNode> seenConceptsOutput = new HashSet<TaxonomyNode>();
		for (String outputVal : s.getOutputs()) {
			TaxonomyNode n = taxonomyMap.get(outputVal);
			s.getTaxonomyOutputs().add(n);

			// Also add output to all parent nodes
			Queue<TaxonomyNode> queue = new LinkedList<TaxonomyNode>();
			queue.add( n );

			while (!queue.isEmpty()) {
			    TaxonomyNode current = queue.poll();
		        seenConceptsOutput.add( current );
		        current.servicesWithOutput.add(s);
		        for (TaxonomyNode parent : current.parents) {
		            if (!seenConceptsOutput.contains( parent )) {
		                queue.add(parent);
		                seenConceptsOutput.add(parent);
		            }
		        }
			}
		}
	}

	private void populateInputs(ServiceNode s) {
		// Populate inputs
		Set<TaxonomyNode> seenConceptsInput = new HashSet<TaxonomyNode>();
		for (String inputVal : s.getInputs()) {
			TaxonomyNode n = taxonomyMap.get(inputVal);
			//n.servicesWithInput.add(s);

			// Also add input to all children nodes
			Queue<TaxonomyNode> queue = new LinkedList<TaxonomyNode>();
			queue.add( n );

			while(!queue.isEmpty()) {
				TaxonomyNode current = queue.poll();
				seenConceptsInput.add( current );

			    Set<String> inputs = current.servicesWithInput.get(s);
			    if (inputs == null) {
			    	inputs = new HashSet<String>();
			    	inputs.add(inputVal);
			    	current.servicesWithInput.put(s, inputs);
			    }
			    else {
			    	inputs.add(inputVal);
			    }

			    for (TaxonomyNode child : current.children) {
			        if (!seenConceptsInput.contains( child )) {
			            queue.add(child);
			            seenConceptsInput.add( child );
			        }
			    }
			}
		}
	}

	/**
	 * Parses the WSC Web service file with the given name, creating Web
	 * services based on this information and saving them to the service map.
	 *
	 * @param fileName
	 */
	private void parseWSCServiceFile(String fileName) {
      Set<String> inputs = new HashSet<String>();
      Set<String> outputs = new HashSet<String>();
      double[] qos = new double[4];

      try {
      	File fXmlFile = new File(fileName);
      	DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
      	DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
      	Document doc = dBuilder.parse(fXmlFile);

      	NodeList nList = doc.getElementsByTagName("service");

      	for (int i = 0; i < nList.getLength(); i++) {
      		org.w3c.dom.Node nNode = nList.item(i);
      		Element eElement = (Element) nNode;

      		String name = eElement.getAttribute("name");
  		    qos[TIME] = Double.valueOf(eElement.getAttribute("Res"));
  		    qos[COST] = Double.valueOf(eElement.getAttribute("Pri"));
  		    qos[AVAILABILITY] = Double.valueOf(eElement.getAttribute("Ava"));
  		    qos[RELIABILITY] = Double.valueOf(eElement.getAttribute("Rel"));


			// Get inputs
			org.w3c.dom.Node inputNode = eElement.getElementsByTagName("inputs").item(0);
			NodeList inputNodes = ((Element)inputNode).getElementsByTagName("instance");
			for (int j = 0; j < inputNodes.getLength(); j++) {
				org.w3c.dom.Node in = inputNodes.item(j);
				Element e = (Element) in;
				inputs.add(e.getAttribute("name"));
			}

			// Get outputs
			org.w3c.dom.Node outputNode = eElement.getElementsByTagName("outputs").item(0);
			NodeList outputNodes = ((Element)outputNode).getElementsByTagName("instance");
			for (int j = 0; j < outputNodes.getLength(); j++) {
				org.w3c.dom.Node out = outputNodes.item(j);
				Element e = (Element) out;
				outputs.add(e.getAttribute("name"));
			}

              ServiceNode ws = new ServiceNode(name, qos, inputs, outputs);
              serviceMap.put(name, ws);
              inputs = new HashSet<String>();
              outputs = new HashSet<String>();
              qos = new double[4];
      	}
      }
      catch(IOException ioe) {
          System.out.println("Service file parsing failed...");
      }
      catch (ParserConfigurationException e) {
          System.out.println("Service file parsing failed...");
		}
      catch (SAXException e) {
          System.out.println("Service file parsing failed...");
		}
  }

	/**
	 * Parses the WSC task file with the given name, extracting input and
	 * output values to be used as the composition task.
	 *
	 * @param fileName
	 */
	private void parseWSCTaskFile(String fileName) {
		try {
	    	File fXmlFile = new File(fileName);
	    	DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
	    	DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
	    	Document doc = dBuilder.parse(fXmlFile);

	    	org.w3c.dom.Node provided = doc.getElementsByTagName("provided").item(0);
	    	NodeList providedList = ((Element) provided).getElementsByTagName("instance");
	    	taskInput = new HashSet<String>();
	    	for (int i = 0; i < providedList.getLength(); i++) {
				org.w3c.dom.Node item = providedList.item(i);
				Element e = (Element) item;
				taskInput.add(e.getAttribute("name"));
	    	}

	    	org.w3c.dom.Node wanted = doc.getElementsByTagName("wanted").item(0);
	    	NodeList wantedList = ((Element) wanted).getElementsByTagName("instance");
	    	taskOutput = new HashSet<String>();
	    	for (int i = 0; i < wantedList.getLength(); i++) {
				org.w3c.dom.Node item = wantedList.item(i);
				Element e = (Element) item;
				taskOutput.add(e.getAttribute("name"));
	    	}
		}
		catch (ParserConfigurationException e) {
          System.out.println("Task file parsing failed...");
          e.printStackTrace();
		}
		catch (SAXException e) {
          System.out.println("Task file parsing failed...");
          e.printStackTrace();
		}
		catch (IOException e) {
          System.out.println("Task file parsing failed...");
          e.printStackTrace();
		}
	}

	/**
	 * Parses the WSC taxonomy file with the given name, building a
	 * tree-like structure.
	 *
	 * @param fileName
	 */
	private void parseWSCTaxonomyFile(String fileName) {
		try {
	    	File fXmlFile = new File(fileName);
	    	DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
	    	DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
	    	Document doc = dBuilder.parse(fXmlFile);
	    	NodeList taxonomyRoots = doc.getChildNodes();

	    	processTaxonomyChildren(null, taxonomyRoots);
		}

		catch (ParserConfigurationException e) {
          System.err.println("Taxonomy file parsing failed...");
		}
		catch (SAXException e) {
          System.err.println("Taxonomy file parsing failed...");
		}
		catch (IOException e) {
          System.err.println("Taxonomy file parsing failed...");
		}
	}

	/**
	 * Recursive function for recreating taxonomy structure from file.
	 *
	 * @param parent - Nodes' parent
	 * @param nodes
	 */
	private void processTaxonomyChildren(TaxonomyNode parent, NodeList nodes) {
		if (nodes != null && nodes.getLength() != 0) {
			for (int i = 0; i < nodes.getLength(); i++) {
				org.w3c.dom.Node ch = nodes.item(i);

				if (!(ch instanceof Text)) {
					Element currNode = (Element) nodes.item(i);
					String value = currNode.getAttribute("name");
					TaxonomyNode taxNode = taxonomyMap.get( value );
					if (taxNode == null) {
					    taxNode = new TaxonomyNode(value);
					    taxonomyMap.put( value, taxNode );
					}
					if (parent != null) {
					    taxNode.parents.add(parent);
						parent.children.add(taxNode);
					}

					NodeList children = currNode.getChildNodes();
					processTaxonomyChildren(taxNode, children);
				}
			}
		}
	}

	private String edgesToDot(Map<String,Edge> edges) {
		StringBuilder builder = new StringBuilder();
		builder.append("digraph g{\n");
		for (Edge e : edges.values()) {
			builder.append(e.toString() + "\n");
		}
		builder.append("}");
		return builder.toString();
	}
}