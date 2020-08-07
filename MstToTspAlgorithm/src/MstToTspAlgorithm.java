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
    	//long exectime;
    	
    	TSPObject(String order,double length){ //,long exectime){
    		this.order = order;
    		this.length = length;
    		//this.exectime = exectime;
    		
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
		/*
		 * public long getExectime() { return exectime;
		 * 
		 * } public void setExectime(long exectime) { this.exectime = exectime;
		 * 
		 * }
		 */
    	
    	
    	public String toString() {
    		return "order: "+order+" length: "+Double.toString(length); //+" execution time: "+Double.toString(exectime);
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
		String inputfile = "distanceMatrix.txt";
		String outputfile = "result.txt";
		
		//randomized: 0 basic version, 1 randomized tree creation
		int randomized = 1;
		int repeats = 100;
		
		try {
			inputfile = args[0];
			randomized = Integer.parseInt(args[1]);
		}
		catch(ArrayIndexOutOfBoundsException e){
    		System.out.println("Not all possible parameters given");
    		
    	}
		
		ArrayList<ArrayList<Double>> distancematrix = new ArrayList<ArrayList<Double>>();
		distancematrix = readfile(inputfile);
		
		if(randomized == 0) {
			//create complete graph
			Pair<ArrayList<Node>, ArrayList<Edge>> vals = createGraph(distancematrix);
			ArrayList<Node> nodes = vals.getFirstValue();
			//ArrayList<Edge> alledges = vals.getEdgelist();
			
			//create minimum spanning tree
			ArrayList<Node> mst = prim(nodes);
			
			//create tsp route
			ArrayList<Node> tsp = mstToTsp(mst);//, alledges, mstedges);
			Pair<String, Double> tspvals = calculateRouteLength(tsp);
			System.out.println("Order: " + tspvals.getFirstValue());
			System.out.println("Length: " + tspvals.getSecondValue());
			
			//write result to txt file
			writeTsp(tspvals.getSecondValue(), tspvals.getFirstValue(), outputfile);
		}
		if(randomized !=0) {
			//randomized version
			ArrayList<TSPObject> list = new ArrayList<TSPObject>();
			
			//create graph
			Pair<ArrayList<Node>, ArrayList<Edge>>vals2 = createGraph(distancematrix);
			ArrayList<Node>nodes2 = vals2.getFirstValue();
			//ArrayList<Edge>alledges2 = vals2.getEdgelist();
			//default 100 repeats
			for (int i = 0; i < repeats; i++) {
				
				//create spanning tree
				ArrayList<Node> st = randomprim(nodes2);
				
				//create tsp route
				ArrayList<Node> tsp2 = mstToTsp(st);//, alledges2, stedges);
				Pair<String, Double> tsp2vals = calculateRouteLength(tsp2);
				//list.add(tsp2vals);
				TSPObject tspobj = new TSPObject(tsp2vals.getFirstValue(), tsp2vals.getSecondValue());
				list.add(tspobj);
				
			}
			Collections.sort(list);
			//System.out.println(list);
			
			TSPObject first = list.get(0);
			System.out.println("Random ver Order: " + first.getOrder());
			System.out.println("Random ver Length: " + first.getLength());
			//write result to txt file
			writeTsp(first.getLength(), first.getOrder(), outputfile);
		}
		

	}

	/**
	 * calculates the length of the tsp and order of the nodes
	 * 
	 * @param tsp list of nodes
	 * @return Pair that has order of the nodes visited and length of the tsp
	 */
	private static Pair<String, Double> calculateRouteLength(ArrayList<Node> tsp) {
		Double routelength = 0.0;
		String order = "";
		ArrayList<Node> processednodes = new ArrayList<Node>();
		int[] node_degree2 = new int[tsp.size()];
		int startnode = 0;
		int nextnode = 0;
		int previousnode = 0;
		
		//calculates node degree for each node (how many edges are connected for each node)
		for (int i = 0; i < tsp.size(); i++) {
			int numofedges = 0;
			for (int j = 0; j < tsp.get(i).edges.size(); j++) {
				if (tsp.get(i).edges.get(j).visited == true) {
					numofedges++;
				}
			}
			node_degree2[i] = numofedges;
		}
		//search leaf node
		for (int i = 0; i < node_degree2.length; i++) {
			if (node_degree2[i] == 1) {
				startnode = i;
				break;
			}
		}
		//find visited edge for the leafnode, add node to order and processed nodes, set leafnode as previous
		for (int i = 0; i < tsp.get(startnode).getEdges().size(); i++) {
			if (tsp.get(startnode).getEdges().get(i).visited == true) {
				routelength = routelength + tsp.get(startnode).getEdges().get(i).length;
				nextnode = tsp.get(startnode).getEdges().get(i).to;
				processednodes.add(tsp.get(startnode));
				previousnode = startnode;
				order = order + previousnode;
			}
		}
		//go through the tsp route to calculate length and order
		while (processednodes.size() < tsp.size()) {
			for (int i = 0; i < tsp.get(nextnode).getEdges().size(); i++) {
				if (tsp.get(nextnode).getEdges().get(i).visited == true
						&& tsp.get(nextnode).getEdges().get(i).to != previousnode) {
					routelength = routelength + tsp.get(nextnode).getEdges().get(i).length;
					processednodes.add(tsp.get(nextnode));
					previousnode = nextnode;
					nextnode = tsp.get(nextnode).getEdges().get(i).to;
					order = order + "-" + previousnode;
				}

			}
			if (node_degree2[nextnode] < 2) {
				processednodes.add(tsp.get(nextnode));
				order = order + "-" + nextnode;
			}
		}
		Pair<String, Double> values = new Pair1<String, Double>(order, routelength);
		return values;
	}

	/**
	 * mst to tsp algorithm, converts minimum spanning tree to open loop travelling
	 * salesman route
	 * 
	 * @param mst      list of nodes
	 * @param edges    list of all edges
	 * @param mstedges list of edges that belong to mst
	 * @return mst list of nodes, with updated edges
	 */
	private static ArrayList<Node> mstToTsp(ArrayList<Node> mst){//, ArrayList<Edge> edges, ArrayList<Edge> mstedges) {
		// ArrayList<Integer> nodeedgecount = new ArrayList<Integer>();
		int[] node_degree = new int[mst.size()];
		// create array with number of edges for each node
		for (int i = 0; i < mst.size(); i++) {
			int numofedges = 0;
			for (int j = 0; j < mst.get(i).edges.size(); j++) {
				if (mst.get(i).edges.get(j).visited == true) {
					numofedges++;
				}
			}
			node_degree[i] = numofedges;

		}
		/*
		 * for (int i = 0; i < mst.size(); i++) { System.out.println("Node: " +
		 * mst.get(i).getNumber() + " Nodedegree: " + node_degree[i]); }
		 */

		int edgetoremove = -1;
		double maxlength = 0.0;
		// going through tree nodes one by one
		for (int i = 0; i < mst.size(); i++) {
			// while node has more edges than 2 (is a knot)
			while (node_degree[i] > 2) {

				// remove edge
				// go through all edges that are connected to node
				for (int j = 0; j < mst.get(i).edges.size(); j++) {
					// if current edge belongs to mst and is longer than previously chosen edge,
					// choose it
					if ((mst.get(i).edges.get(j).visited == true) && (mst.get(i).edges.get(j).length > maxlength)) {
						maxlength = mst.get(i).edges.get(j).length;
						edgetoremove = j;
					}
				}
				maxlength = 0.0;
				
				// set removed edge to not visited for current node
				mst.get(i).edges.get(edgetoremove).visited = false;

				// update edgecount for the node
				node_degree[i] = node_degree[i] - 1;

				// Choose the node that was on the other end of the deleted edge
				int otherside = mst.get(i).edges.get(edgetoremove).to;
				
				// Go through all the edges for that node
				for (int j = 0; j < mst.get(otherside).edges.size(); j++) {
					// remove the same edge for that node too
					if (mst.get(otherside).edges.get(j).to == mst.get(i).getNumber()) {
						mst.get(otherside).edges.get(j).visited = false;
					}

				}
				// update nodecount for that node too
				node_degree[otherside] = node_degree[otherside] - 1;

				// edge succesfully removed from mst

				// Create set for the subgraph that was disconnected
				HashSet<Integer> set = new HashSet<>();
				set.add(otherside);
				for (int j = 0; j < mst.get(otherside).edges.size(); j++) {
					if (mst.get(otherside).edges.get(j).visited == true) {
						set.add(mst.get(otherside).edges.get(j).to);
						Set<Integer> subset = addSet(mst, mst.get(otherside).edges.get(j).to, set);
						set.addAll(subset);
					}
				}
				
				// add edge to otherside node
				int tonode = 0;
				// go through all edges of the node
				for (int j = 0; j < mst.get(otherside).edges.size(); j++) {
					// add edge to between this node and another leafnode // check also that
					// leafnode is not in the same subgraph
					if (mst.get(otherside).edges.get(j).visited == false
							&& node_degree[mst.get(otherside).edges.get(j).to] < 2
							&& !set.contains(mst.get(otherside).edges.get(j).to)) {

						mst.get(otherside).edges.get(j).visited = true;

						// update edgecount for this node
						node_degree[otherside] = node_degree[otherside] + 1;

						// the leaf node we are connecting the edge
						tonode = mst.get(otherside).edges.get(j).to;
						break;

					}
				}
				// go through all the edges of the leafnode we were connecting the edge
				for (int j = 0; j < mst.get(tonode).edges.size(); j++) {
					// Choose the edge for this node where we came
					if (mst.get(tonode).edges.get(j).to == otherside) {
						mst.get(tonode).edges.get(j).visited = true;
						node_degree[tonode] = node_degree[tonode] + 1;
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
				HashSet<Integer> subset = addSet(mst, mst.get(to).edges.get(j).to, set);
				set.addAll(subset);
			}
		}

		return set;
	}

	/**
	 * prims algorithm
	 * 
	 * @param nodes list of nodes
	 * @param edges list of all edges
	 * @return values pair that contains nodes and edges that belong to mst
	 */
	private static ArrayList<Node> prim(ArrayList<Node> nodes){//, ArrayList<Edge> edges) {
		ArrayList<Node> mst = new ArrayList<Node>();
		ArrayList<Edge> mstedges = new ArrayList<Edge>();

		int startnode = (int) (Math.random() * nodes.size());
		nodes.get(startnode).chosen = true;
		mst.add(nodes.get(startnode));

		Double shortest = null;
		Edge edgetoadd = null;
		int edgeindex = 0;

		//while not every node is connected to the tree
		while (mst.size() < nodes.size()) {
			//go through every node that belongs to the tree
			for (int i = 0; i < mst.size(); i++) {
				//go through all of the edges for those nodes
				for (int j = 0; j < mst.get(i).edges.size(); j++) {
					//if the edge is not visited and the destination node is not already chosen, edge is one possibility
					if (mst.get(i).edges.get(j).visited == false && nodes.get(mst.get(i).edges.get(j).to).chosen == false) {
						
						Edge possibleedge = mst.get(i).edges.get(j);
						if (edgetoadd == null) {
							shortest = possibleedge.getLength();
							edgetoadd = possibleedge;
							edgeindex = j;
						}
						//if current possible is shorter than previous, update
						if (possibleedge.getLength() < shortest) {
							shortest = possibleedge.getLength();
							edgetoadd = possibleedge;
							edgeindex = j;
						}
					}
				}
			}
			nodes.get(edgetoadd.to).chosen = true;
			//update also the same edge that belongs to the destination node
			for (int i = 0; i < nodes.get(edgetoadd.to).edges.size(); i++) {
				if (nodes.get(edgetoadd.to).edges.get(i).to == edgetoadd.from) {
					nodes.get(edgetoadd.to).edges.get(i).visited = true;
				}
			}
			nodes.get(edgetoadd.from).edges.get(edgeindex).visited = true;
			
			mstedges.add(edgetoadd);

			mst.add(nodes.get(edgetoadd.to));
			edgetoadd = null;
			
		}
		Collections.sort(mst);
		
		return mst;
		//Pair<ArrayList<Node>, ArrayList<Edge>> values = new Pair1<ArrayList<Node>, ArrayList<Edge>>(mst, mstedges);
		//return values;
	}

	/**
	 * randomized tree creation, creates spannigtree
	 * prim with randomization, chooses edge to add randomly from 3 shortest eligble edges
	 * 
	 * @param nodes list of nodes
	 * @param edges list of all edges
	 * @return values pair that contains nodes and edges that belong to spanningtree
	 */
	private static ArrayList<Node> randomprim(ArrayList<Node> nodes){//, ArrayList<Edge> edges) {
		
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
		int numofshortest = 3;

		int startnode = (int) (Math.random() * nodes.size());
		nodes.get(startnode).chosen = true;
		st.add(nodes.get(startnode));
		
		Edge edgetoadd = null;
		
		//while not every node is connected to the tree
		while (st.size() < nodes.size()) {
			//go through every node that belongs to the tree
			for (int i = 0; i < st.size(); i++) {
				//go through all of the edges for those nodes
				for (int j = 0; j < st.get(i).edges.size(); j++) {
					//if the edge is not visited and the destination node is not already chosen, edge is one possibility
					if (st.get(i).edges.get(j).visited == false && nodes.get(st.get(i).edges.get(j).to).chosen == false) {
						Edge possibleedge = st.get(i).edges.get(j);
						possible.add(possibleedge);
						
					}
				}
			}
			//sort possibilities
			Collections.sort(possible);
			int index = 0;

			//choose randomly the edge from shortest 3 or if less possibilities then from all of them (2 or 1)
			if (possible.size() > numofshortest) {
				index = randomizer.nextInt(numofshortest);
			}
			else {
				index = randomizer.nextInt(possible.size());
			} 
			edgetoadd = possible.get(index);

			nodes.get(edgetoadd.to).chosen = true;
			for (int i = 0; i < nodes.get(edgetoadd.to).edges.size(); i++) {
				if (nodes.get(edgetoadd.to).edges.get(i).to == edgetoadd.from) {
					nodes.get(edgetoadd.to).edges.get(i).visited = true;
					
				}
			}
			for (int i = 0; i < nodes.get(edgetoadd.from).edges.size(); i++) {
				if (nodes.get(edgetoadd.from).edges.get(i).to == edgetoadd.to) {
					nodes.get(edgetoadd.from).edges.get(i).visited = true;
				}
			}
			//stedges.add(edgetoadd);

			st.add(nodes.get(edgetoadd.to));
			edgetoadd = null;
			possible.clear();
		}

		Collections.sort(st);
		return st;
		//Pair<ArrayList<Node>, ArrayList<Edge>> values = new Pair1<ArrayList<Node>, ArrayList<Edge>>(st, stedges);
		//return values;

	}

	/**
	 * writes result to specified file
	 * @param total_length length of the tsp
	 * @param order        order of the nodes
	 * @param outputfile   file where result is written
	 * @throws IOException
	 */
	private static void writeTsp(Double total_length, String order, String outputfile) throws IOException {
		PrintWriter writer = new PrintWriter(new FileWriter(outputfile));
		writer.println(order);
		writer.println();
		writer.print("&" + total_length);
		writer.close();
	}
	
	/**
	 * reads file from specified path (distancematrix in txt file)
	 * @param inputfile filepath
	 * @return matrix distancematrix
	 * @throws IOException
	 */
	private static ArrayList<ArrayList<Double>> readfile(String inputfile) throws IOException {
		FileReader in = new FileReader(inputfile);
		LineNumberReader lineNumberReader = new LineNumberReader(new BufferedReader(in));
		
		ArrayList<ArrayList<Double>> matrix = new ArrayList<ArrayList<Double>>();
		
		String currentline;
		while((currentline = lineNumberReader.readLine()) != null) {
        	if(currentline.contains("EOF")) {
        		break;
        	}
        	String[] nodes = currentline.split(",");
        	ArrayList<Double> row = new ArrayList<Double>();
        	for(int i = 0; i < nodes.length;i++) {
        		row.add(Double.parseDouble(nodes[i]));
        	}
        	matrix.add(row);
		
		}
		return matrix;
	}
	
	

	// Pair interface
	public interface Pair<A, B> { 
		public A getFirstValue();
		public B getSecondValue(); 
	}
	// Pair class
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
	 * @param distancematrix matrix that has distances between every node
	 * @return values pair that has nodes and edges that belong to the graph
	 */
	private static Pair<ArrayList<Node>, ArrayList<Edge>> createGraph(ArrayList<ArrayList<Double>> distancematrix) {
		ArrayList<Node> nodelist = new ArrayList<Node>();
		ArrayList<Edge> edgelist = new ArrayList<Edge>();
		ArrayList<String> edgeswaplist = new ArrayList<String>();
		int numofnodes = distancematrix.size();
		for (int i = 0; i < numofnodes; i++) {
			nodelist.add(new Node(i));
		}
		for (int i = 0; i < distancematrix.size(); i++) {
			for (int j = 0; j < distancematrix.size(); j++) {
				if (i != j) {
					Edge edge = new Edge(i, j, distancematrix.get(i).get(j));
					nodelist.get(i).addEdge(edge);
					String edgeswap = edge.to + "," + edge.from;
					String edgeswap2 = edge.from + "," + edge.to;
					if (!edgeswaplist.contains(edgeswap) && !edgeswaplist.contains(edgeswap2)) {
						edgeswaplist.add(edgeswap);
						edgeswaplist.add(edgeswap2);
						edgelist.add(edge);
					}
				}
			}
		}
		for (int i = 0; i < nodelist.size(); i++) {
			Collections.sort(nodelist.get(i).edges);
		}
		Collections.sort(edgelist);
		Pair<ArrayList<Node>, ArrayList<Edge>> values = new Pair1<ArrayList<Node>, ArrayList<Edge>>(nodelist, edgelist);
		return values;
	}

}
