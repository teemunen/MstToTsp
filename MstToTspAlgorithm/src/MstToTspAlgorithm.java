import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * MstToTsp algorithm
 * Solves openloop travelling salesman problem by modifying minimum spanningtree
 * Also includes randomized version that creates randomized spanning tree and modifies that for tsp route
 * @author Teemu Nenonen
 *
 */
public class MstToTspAlgorithm {
	/**
	 * Class used to store node information
	 */
	private static class Node implements Comparable<Node> {
		int number;
		boolean chosen;
		ArrayList<Edge> edges;
		/**
		 * Default constructor of a node, chosen is false and edgelist for the node is emty
		 * @param number number of the node
		 */
		Node(int number) {
			this.number = number;
			this.chosen = false;
			this.edges = new ArrayList<Edge>();
		}
		public int getNumber() {
			return number;
		}
		public void setNumber(int number) {
			this.number = number;
		}
		public boolean isChosen() {
			return chosen;
		}
		public void setChosen(boolean chosen) {
			this.chosen = chosen;
		}
		public ArrayList<Edge> getEdges() {
			return edges;
		}
		public void setEdges(ArrayList<Edge> edges) {
			this.edges = edges;
		}
		public void addEdge(Edge edge) {
			this.edges.add(edge);
		}
		public String toString() {
			return "Node: " + number + " chosen: " + chosen + " edgelist: " + edges;
		}
		/**
		 * Override of the compareTo method, used to sort nodes by number
		 * @param node node to compare
		 * @return 0 if their number is the same, 1 if current node number is more than compared one, and -1 otherwise
		 */
		@Override
		public int compareTo(Node node) {
			int value = 0;
			int num = node.getNumber();
			if (this.number > num) {
				value = 1;
			} else if (this.number < num) {
				value = -1;
			} else if (this.number == num) {
				value = 0;
			}
			return value;
		}
	}
	/**
	 * Class used to store information of edge
	 */
	private static class Edge implements Comparable<Edge> {
		int from;
		int to;
		double length;
		boolean visited;
		/**
		 * Default constructor of an edge
		 * @param from edge startpoint
		 * @param to edge endpoint
		 * @param length length of the edge
		 */
		Edge(int from, int to, double length) {
			this.from = from;
			this.to = to;
			this.length = length;
			this.visited = false;
		}
		public double getLength() {
			return length;
		}
		public void setLength(double length) {
			this.length = length;
		}
		/**
		 * Override of the compareTo method, used to sort edges
		 * @param edge edge to compare
		 * @return 0 if their lengths are the same, 1 if current edge is longer than compared one, and -1 otherwise
		 */
		@Override
		public int compareTo(Edge edge) {
			int value = 0;
			Double compareLength = edge.getLength();
			if (this.length > compareLength) {
				value = 1;
			} else if (this.length < compareLength) {
				value = -1;
			} else if (this.length == compareLength) {
				value = 0;
			}
			return value;
		}
		public String toString() {
			return "From: " + from + " To: " + to + " length: " + length + " visited: " + visited;
		}
	}
	/**
     * Class used to store TSP information
     */
    private static class TSPObject implements Comparable<TSPObject>{
    	String order;
    	double length;
    	
    	TSPObject(String order,double length){ 
    		this.order = order;
    		this.length = length;
    	}
    	public String getOrder() {
			return order;
    	}
    	public void setOrder(String order) {
    		this.order = order;
    	}
    	public Double getLength() {
			return length;
    	}
    	public void setLength(double length) {
    		this.length = length;
    	}
    	public String toString() {
    		return "order: "+order+" length: "+Double.toString(length);
    	}
    	/**
    	 * Override of the compareTo method, used to sort tsp routes by length 
         * @param TSPObject tspobj to compare
         * 
    	 */
		@Override
		public int compareTo(TSPObject tspobj) {
			int value = 0;
			Double compareLength = tspobj.getLength();
			if(this.length>compareLength) {
				value = 1;	
			}
			else if(this.length<compareLength){
				value = -1;
			}
			else if(this.length==compareLength){
				value = 0;
			}
			return value;
		}
    }
	public static void main(String[] args) throws IOException{
		//Default filepaths
		String inputFile = "distanceMatrix.txt";
		String outputFile = "result.txt";
		
		//randomized: 0 basic version, 1 randomized tree creation
		int randomized = 1;
		int repeats = 100;
		
		//If given filepaths as parameter then uses those
		try {
			inputFile = args[0];
			randomized = Integer.parseInt(args[1]);
		}
		catch(ArrayIndexOutOfBoundsException e){
    		System.out.println("Not all possible parameters given");
    	}
		ArrayList<ArrayList<Double>> distanceMatrix = new ArrayList<ArrayList<Double>>();
		distanceMatrix = readFile(inputFile);
		
		if(randomized == 0) {
			//create complete graph
			Pair<ArrayList<Node>, ArrayList<Edge>> vals = createGraph(distanceMatrix);
			ArrayList<Node> nodes = vals.getFirstValue();
			//ArrayList<Edge> alledges = vals.getEdgelist();
			
			//create minimum spanning tree
			ArrayList<Node> mst = prim(nodes);
			
			//create tsp route
			ArrayList<Node> tsp = mstToTsp(mst);//, alledges, mstedges);
			Pair<String, Double> tspVals = calculateRouteLength(tsp);
			System.out.println("Order: " + tspVals.getFirstValue());
			System.out.println("Length: " + tspVals.getSecondValue());
			
			//write result to txt file
			writeTspResultToFile(tspVals.getSecondValue(), tspVals.getFirstValue(), outputFile);
		}
		if(randomized !=0) {
			//randomized version
			ArrayList<TSPObject> list = new ArrayList<TSPObject>();
			
			//create graph
			Pair<ArrayList<Node>, ArrayList<Edge>>vals2 = createGraph(distanceMatrix);
			ArrayList<Node>nodes2 = vals2.getFirstValue();
			//ArrayList<Edge>alledges2 = vals2.getEdgelist();
			//default 100 repeats
			for (int i = 0; i < repeats; i++) {
				//create spanning tree
				ArrayList<Node> st = randomPrim(nodes2);
				
				//create tsp route
				ArrayList<Node> tsp2 = mstToTsp(st);//, alledges2, stedges);
				Pair<String, Double> tsp2Vals = calculateRouteLength(tsp2);
				//list.add(tsp2vals);
				TSPObject tspobj = new TSPObject(tsp2Vals.getFirstValue(), tsp2Vals.getSecondValue());
				list.add(tspobj);
			}
			Collections.sort(list);
			//System.out.println(list);
			
			TSPObject first = list.get(0);
			System.out.println("Random ver Order: " + first.getOrder());
			System.out.println("Random ver Length: " + first.getLength());
			//write result to txt file
			writeTspResultToFile(first.getLength(), first.getOrder(), outputFile);
		}
	}
	/**
	 * calculates the length of the tsp and order of the nodes
	 * 
	 * @param tsp list of nodes
	 * @return Pair that has order of the nodes visited and length of the tsp
	 */
	private static Pair<String, Double> calculateRouteLength(ArrayList<Node> tsp) {
		Double routeLength = 0.0;
		String order = "";
		ArrayList<Node> processedNodes = new ArrayList<Node>();
		int[] nodeDegree2 = new int[tsp.size()];
		int startNode = 0;
		int nextNode = 0;
		int previousNode = 0;
		
		//calculates node degree for each node (how many edges are connected for each node)
		for (int i = 0; i < tsp.size(); i++) {
			int numofedges = 0;
			for (int j = 0; j < tsp.get(i).edges.size(); j++) {
				if (tsp.get(i).edges.get(j).visited == true) {
					numofedges++;
				}
			}
			nodeDegree2[i] = numofedges;
		}
		//search leaf node
		for (int i = 0; i < nodeDegree2.length; i++) {
			if (nodeDegree2[i] == 1) {
				startNode = i;
				break;
			}
		}
		//find visited edge for the leafnode, add node to order and processed nodes, set leafnode as previous
		for (int i = 0; i < tsp.get(startNode).getEdges().size(); i++) {
			if (tsp.get(startNode).getEdges().get(i).visited == true) {
				routeLength = routeLength + tsp.get(startNode).getEdges().get(i).length;
				nextNode = tsp.get(startNode).getEdges().get(i).to;
				processedNodes.add(tsp.get(startNode));
				previousNode = startNode;
				order = order + previousNode;
			}
		}
		//go through the tsp route to calculate length and order
		while (processedNodes.size() < tsp.size()) {
			for (int i = 0; i < tsp.get(nextNode).getEdges().size(); i++) {
				if (tsp.get(nextNode).getEdges().get(i).visited == true
						&& tsp.get(nextNode).getEdges().get(i).to != previousNode) {
					routeLength = routeLength + tsp.get(nextNode).getEdges().get(i).length;
					processedNodes.add(tsp.get(nextNode));
					previousNode = nextNode;
					nextNode = tsp.get(nextNode).getEdges().get(i).to;
					order = order + "-" + previousNode;
				}
			}
			if (nodeDegree2[nextNode] < 2) {
				processedNodes.add(tsp.get(nextNode));
				order = order + "-" + nextNode;
			}
		}
		Pair<String, Double> values = new Pair1<String, Double>(order, routeLength);
		return values;
	}
	/**
	 * mst to tsp algorithm, converts minimum spanning tree to open loop travelling
	 * salesman route
	 * 
	 * @param mst      list of nodes
	 * @return mst list of nodes, with updated edges
	 */
	private static ArrayList<Node> mstToTsp(ArrayList<Node> mst){
		// ArrayList<Integer> nodeedgecount = new ArrayList<Integer>();
		int[] nodeDegree = new int[mst.size()];
		// create array with number of edges for each node
		for (int i = 0; i < mst.size(); i++) {
			int numOfEdges = 0;
			for (int j = 0; j < mst.get(i).edges.size(); j++) {
				if (mst.get(i).edges.get(j).visited == true) {
					numOfEdges++;
				}
			}
			nodeDegree[i] = numOfEdges;
		}
		/*
		 * for (int i = 0; i < mst.size(); i++) { System.out.println("Node: " +
		 * mst.get(i).getNumber() + " Nodedegree: " + node_degree[i]); }
		 */
		int edgeToRemove = -1;
		double maxLength = 0.0;
		// going through tree nodes one by one
		for (int i = 0; i < mst.size(); i++) {
			// while node has more edges than 2 (is a knot)
			while (nodeDegree[i] > 2) {
				// remove edge
				// go through all edges that are connected to node
				for (int j = 0; j < mst.get(i).edges.size(); j++) {
					// if current edge belongs to mst and is longer than previously chosen edge,
					// choose it
					if ((mst.get(i).edges.get(j).visited == true) && (mst.get(i).edges.get(j).length > maxLength)) {
						maxLength = mst.get(i).edges.get(j).length;
						edgeToRemove = j;
					}
				}
				maxLength = 0.0;
				
				// set removed edge to not visited for current node
				mst.get(i).edges.get(edgeToRemove).visited = false;

				// update edgecount for the node
				nodeDegree[i] = nodeDegree[i] - 1;

				// Choose the node that was on the other end of the deleted edge
				int otherSide = mst.get(i).edges.get(edgeToRemove).to;
				
				// Go through all the edges for that node
				for (int j = 0; j < mst.get(otherSide).edges.size(); j++) {
					// remove the same edge for that node too
					if (mst.get(otherSide).edges.get(j).to == mst.get(i).getNumber()) {
						mst.get(otherSide).edges.get(j).visited = false;
					}
				}
				// update nodecount for that node too
				nodeDegree[otherSide] = nodeDegree[otherSide] - 1;

				// edge succesfully removed from mst

				// Create set for the subgraph that was disconnected
				HashSet<Integer> set = new HashSet<>();
				set.add(otherSide);
				for (int j = 0; j < mst.get(otherSide).edges.size(); j++) {
					if (mst.get(otherSide).edges.get(j).visited == true) {
						set.add(mst.get(otherSide).edges.get(j).to);
						Set<Integer> subSet = addSet(mst, mst.get(otherSide).edges.get(j).to, set);
						set.addAll(subSet);
					}
				}
				// add edge to otherside node
				int toNode = 0;
				// go through all edges of the node
				for (int j = 0; j < mst.get(otherSide).edges.size(); j++) {
					// add edge to between this node and another leafnode // check also that
					// leafnode is not in the same subgraph
					if (mst.get(otherSide).edges.get(j).visited == false
							&& nodeDegree[mst.get(otherSide).edges.get(j).to] < 2
							&& !set.contains(mst.get(otherSide).edges.get(j).to)) {

						mst.get(otherSide).edges.get(j).visited = true;

						// update edgecount for this node
						nodeDegree[otherSide] = nodeDegree[otherSide] + 1;

						// the leaf node we are connecting the edge
						toNode = mst.get(otherSide).edges.get(j).to;
						break;
					}
				}
				// go through all the edges of the leafnode we were connecting the edge
				for (int j = 0; j < mst.get(toNode).edges.size(); j++) {
					// Choose the edge for this node where we came
					if (mst.get(toNode).edges.get(j).to == otherSide) {
						mst.get(toNode).edges.get(j).visited = true;
						nodeDegree[toNode] = nodeDegree[toNode] + 1;
					}
				}
			}
		}
		/*
		 * for (int i = 0; i < mst.size(); i++) { System.out.println("Node: " +
		 * mst.get(i).getNumber() + " Nodedegree: " + node_degree[i]); }
		 */
		return mst;
	}
	/**
	 * 
	 * @param mst list of nodes
	 * @param to number of the node
	 * @param set nodes that are added to set
	 * @return set
	 */
	private static HashSet<Integer> addSet(ArrayList<Node> mst, int to, HashSet<Integer> set) {
		for (int j = 0; j < mst.get(to).edges.size(); j++) {
			if (mst.get(to).edges.get(j).visited == true && !set.contains(mst.get(to).edges.get(j).to)) {
				set.add(mst.get(to).edges.get(j).to);
				HashSet<Integer> subSet = addSet(mst, mst.get(to).edges.get(j).to, set);
				set.addAll(subSet);
			}
		}
		return set;
	}
	/**
	 * prims algorithm
	 * 
	 * @param nodes list of nodes
	 * @return values pair that contains nodes and edges that belong to mst
	 */
	private static ArrayList<Node> prim(ArrayList<Node> nodes){
		ArrayList<Node> mst = new ArrayList<Node>();
		ArrayList<Edge> mstEdges = new ArrayList<Edge>();

		int startNode = (int) (Math.random() * nodes.size());
		nodes.get(startNode).chosen = true;
		mst.add(nodes.get(startNode));

		Double shortest = null;
		Edge edgeToAdd = null;
		int edgeIndex = 0;

		//while not every node is connected to the tree
		while (mst.size() < nodes.size()) {
			//go through every node that belongs to the tree
			for (int i = 0; i < mst.size(); i++) {
				//go through all of the edges for those nodes
				for (int j = 0; j < mst.get(i).edges.size(); j++) {
					//if the edge is not visited and the destination node is not already chosen, edge is one possibility
					if (mst.get(i).edges.get(j).visited == false && nodes.get(mst.get(i).edges.get(j).to).chosen == false) {
						Edge possibleEdge = mst.get(i).edges.get(j);
						if (edgeToAdd == null) {
							shortest = possibleEdge.getLength();
							edgeToAdd = possibleEdge;
							edgeIndex = j;
						}
						//if current possible is shorter than previous, update
						if (possibleEdge.getLength() < shortest) {
							shortest = possibleEdge.getLength();
							edgeToAdd = possibleEdge;
							edgeIndex = j;
						}
					}
				}
			}
			nodes.get(edgeToAdd.to).chosen = true;
			//update also the same edge that belongs to the destination node
			for (int i = 0; i < nodes.get(edgeToAdd.to).edges.size(); i++) {
				if (nodes.get(edgeToAdd.to).edges.get(i).to == edgeToAdd.from) {
					nodes.get(edgeToAdd.to).edges.get(i).visited = true;
				}
			}
			nodes.get(edgeToAdd.from).edges.get(edgeIndex).visited = true;
			
			mstEdges.add(edgeToAdd);

			mst.add(nodes.get(edgeToAdd.to));
			edgeToAdd = null;
		}
		Collections.sort(mst);
		return mst;
	}
	/**
	 * randomized tree creation, creates spannigtree
	 * prim with randomization, chooses edge to add randomly from 3 shortest eligble edges
	 * 
	 * @param nodes list of nodes
	 * @return values pair that contains nodes and edges that belong to spanningtree
	 */
	private static ArrayList<Node> randomPrim(ArrayList<Node> nodes){
		
		//set every node to not chosen and every edge also to not visited
		for(int j = 0; j < nodes.size(); j++) {
			nodes.get(j).chosen = false;
			for(int k = 0; k < nodes.get(j).edges.size(); k++) {
				nodes.get(j).edges.get(k).visited = false;
			}
		}
		ArrayList<Node> st = new ArrayList<Node>();
		//ArrayList<Edge> stedges = new ArrayList<Edge>();
		ArrayList<Edge> possible = new ArrayList<Edge>();

		Random randomizer = new Random();
		int numOfShortest = 3;

		int startNode = (int) (Math.random() * nodes.size());
		nodes.get(startNode).chosen = true;
		st.add(nodes.get(startNode));
		
		Edge edgeToAdd = null;
		
		//while not every node is connected to the tree
		while (st.size() < nodes.size()) {
			//go through every node that belongs to the tree
			for (int i = 0; i < st.size(); i++) {
				//go through all of the edges for those nodes
				for (int j = 0; j < st.get(i).edges.size(); j++) {
					//if the edge is not visited and the destination node is not already chosen, edge is one possibility
					if (st.get(i).edges.get(j).visited == false && nodes.get(st.get(i).edges.get(j).to).chosen == false) {
						Edge possibleEdge = st.get(i).edges.get(j);
						possible.add(possibleEdge);
					}
				}
			}
			//sort possibilities
			Collections.sort(possible);
			int index = 0;

			//choose randomly the edge from shortest 3 or if less possibilities then from all of them (2 or 1)
			if (possible.size() > numOfShortest) {
				index = randomizer.nextInt(numOfShortest);
			}
			else {
				index = randomizer.nextInt(possible.size());
			} 
			edgeToAdd = possible.get(index);

			nodes.get(edgeToAdd.to).chosen = true;
			for (int i = 0; i < nodes.get(edgeToAdd.to).edges.size(); i++) {
				if (nodes.get(edgeToAdd.to).edges.get(i).to == edgeToAdd.from) {
					nodes.get(edgeToAdd.to).edges.get(i).visited = true;
				}
			}
			for (int i = 0; i < nodes.get(edgeToAdd.from).edges.size(); i++) {
				if (nodes.get(edgeToAdd.from).edges.get(i).to == edgeToAdd.to) {
					nodes.get(edgeToAdd.from).edges.get(i).visited = true;
				}
			}
			//stedges.add(edgetoadd);
			st.add(nodes.get(edgeToAdd.to));
			edgeToAdd = null;
			possible.clear();
		}
		Collections.sort(st);
		return st;
	}
	/**
	 * writes result to specified file
	 * @param totalLength length of the tsp
	 * @param order        order of the nodes
	 * @param outputFile   file where result is written
	 * @throws IOException
	 */
	private static void writeTspResultToFile(Double totalLength, String order, String outputFile) throws IOException {
		PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
		writer.println(order);
		writer.println();
		writer.print("&" + totalLength);
		writer.close();
	}
	/**
	 * reads file from specified path (distancematrix in txt file)
	 * @param inputFile filepath
	 * @return matrix distancematrix
	 * @throws IOException
	 */
	private static ArrayList<ArrayList<Double>> readFile(String inputFile) throws IOException {
		FileReader in = new FileReader(inputFile);
		LineNumberReader lineNumberReader = new LineNumberReader(new BufferedReader(in));
		
		ArrayList<ArrayList<Double>> matrix = new ArrayList<ArrayList<Double>>();
		
		String currentLine;
		while((currentLine = lineNumberReader.readLine()) != null) {
        	if(currentLine.contains("EOF")) {
        		break;
        	}
        	String[] nodes = currentLine.split(",");
        	ArrayList<Double> row = new ArrayList<Double>();
        	for(int i = 0; i < nodes.length;i++) {
        		row.add(Double.parseDouble(nodes[i]));
        	}
        	matrix.add(row);
		}
		lineNumberReader.close();
		return matrix;
	}
	/**
	 *  Pair interface
	 */
	public interface Pair<A, B> { 
		public A getFirstValue();
		public B getSecondValue(); 
	}
	/**
	 *  Pair class
	 */
	public static class Pair1<A, B> implements Pair<A, B> {

		private A first;
		private B second;

		public Pair1(A first, B second) {
			this.first = first;
			this.second = second;
		}
		public A getFirstValue() {
			return first;
		}
		public B getSecondValue() {
			return second;
		}
	}
	/**
	 * Creates complete graph
	 * 
	 * @param distanceMatrix matrix that has distances between every node
	 * @return values pair that has nodes and edges that belong to the graph
	 */
	private static Pair<ArrayList<Node>, ArrayList<Edge>> createGraph(ArrayList<ArrayList<Double>> distanceMatrix) {
		ArrayList<Node> nodeList = new ArrayList<Node>();
		ArrayList<Edge> edgeList = new ArrayList<Edge>();
		ArrayList<String> edgeSwapList = new ArrayList<String>();
		int numOfNodes = distanceMatrix.size();
		for (int i = 0; i < numOfNodes; i++) {
			nodeList.add(new Node(i));
		}
		for (int i = 0; i < distanceMatrix.size(); i++) {
			for (int j = 0; j < distanceMatrix.size(); j++) {
				if (i != j) {
					Edge edge = new Edge(i, j, distanceMatrix.get(i).get(j));
					nodeList.get(i).addEdge(edge);
					String edgeswap = edge.to + "," + edge.from;
					String edgeswap2 = edge.from + "," + edge.to;
					if (!edgeSwapList.contains(edgeswap) && !edgeSwapList.contains(edgeswap2)) {
						edgeSwapList.add(edgeswap);
						edgeSwapList.add(edgeswap2);
						edgeList.add(edge);
					}
				}
			}
		}
		for (int i = 0; i < nodeList.size(); i++) {
			Collections.sort(nodeList.get(i).edges);
		}
		Collections.sort(edgeList);
		Pair<ArrayList<Node>, ArrayList<Edge>> values = new Pair1<ArrayList<Node>, ArrayList<Edge>>(nodeList, edgeList);
		return values;
	}
}
