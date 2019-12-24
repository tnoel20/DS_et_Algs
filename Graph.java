import java.util.NoSuchElementException;
import java.util.List;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.PriorityQueue;
import java.util.ArrayDeque;
import java.util.Collections;
import java.util.Random;
import java.util.LinkedList;
import java.io.IOException;
import java.io.InputStream;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParserFactory;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;


/**
 * A Simple Graph that can be Directed or Undirected
 * 
 * @param <T> The type of the graph's vertices
 * @author The Purple Teletubby
 */
public class Graph<T> 
{
    private static final int UNDEFINED = -1;
    
    private boolean _isDirected;
    
    // The graph's implementation; vertices are mapped to a map containing 
    // both the destination vertex and the edge's weight
    private HashMap<T, HashMap<T, Integer>> _adjacencyMap;
    
    /**
     * Constructs a directed or undirected graph
     * 
     * @param isDirected Indicates whether the graph is directed
     */
    public Graph(boolean isDirected) 
    {
        // Indicates whether or not the graph will be directed
        _isDirected = isDirected;
        _adjacencyMap = new HashMap<>();
    }
    
    
    /**
     * Returns a list of the vertices in a graph
     * 
     * @return A list of vertices
     */
    public List<T> getVertices() 
    {
        // Returns a list containing the graph's vertices
        return new ArrayList<>(_adjacencyMap.keySet());
    }
    
    
    /**
     * Returns the weight of a specified edge
     * 
     * @param source The source vertex of the edge
     * @param destination The destination vertex of the edge
     * @return An integer describing the weight of the edge
     */
    public int getEdgeWeight(T source, T destination) 
    {
        // Weight of UNDEFINED indicates that no such edge exists
        int weight = UNDEFINED;
        
        // If the source and destination are currently connected
        // vertices in the graph
        if (edgeExists(source, destination)) 
        {
            weight = _adjacencyMap.get(source).get(destination);
        }
        
        return weight;
    }
    
    
    /**
     * Gets the edge corresponding to the specified source and
     * destination vertices
     * 
     * @param source The source vertex of the edge
     * @param destination The destination vertex of the edge
     * @return The edge corresponding to the two specified vertices
     */
    public Edge<T> getEdge(T source, T destination) 
    {
        Edge<T> edge;
        
        // If the specifed edge does not exist, the edge is null
        if (!edgeExists(source, destination))
        {
            edge = null;
        }
        else
        {
            edge = new Edge<>(source, destination, 
                    _adjacencyMap.get(source).get(destination));
        }
        
        // Returns an edge corresponding to the specified source and 
        // destination vertices if it exists
        return edge;
    }
    
    
    /**
     * Returns all of the edges in the graph as a list
     * 
     * @return A list of the edges in the graph
     */
    public List<Edge<T>> getEdges() 
    {
        HashSet<T> vertSet = new HashSet<>();
        ArrayList<Edge<T>> edgeList = new ArrayList<>();
        
        for (T sourceVert: _adjacencyMap.keySet())
        {
            // For each adjacent vertex of the sourceVert
            for (T destVert: _adjacencyMap.get(sourceVert).keySet())
            {
                // If the graph is directed, or the destination vertex 
                // has not been a previous source 
                if (_isDirected || !vertSet.contains(destVert))
                {
                    // Add the edge corresponding to the specified source
                    // and destination vertices
                    edgeList.add(getEdge(sourceVert, destVert));
                }
            }

            vertSet.add(sourceVert);
        }
        
        return edgeList;
    }
    
    
    /**
     * Returns whether or not the specified edge exists
     * 
     * @param source The source vertex of the edge
     * @param destination The destination vertex of the edge
     * @return A boolean value that is true if and only if the edge exists
     */
    public boolean edgeExists(T source, T destination) 
    {
        // Returns true if edge exists
        return _adjacencyMap.containsKey(source) && 
                _adjacencyMap.get(source).containsKey(destination);
    }
    
    
    /**
     * Adds a vertex to a graph
     * 
     * @param vertex The vertex to be added to the graph
     * @throws IllegalArgumentException if the vertex already exists
     * in the graph
     */
    public void addVertex(T vertex) throws IllegalArgumentException 
    {
        // If vertex is already one of the graph's vertices
        if (vertex == null || _adjacencyMap.containsKey(vertex)) 
        {
            throw new IllegalArgumentException("ERROR: null vertices and "
                    + "vertices that already exist in graph cannot be added.");
        }
        
        // Adds the new vertex to the graph
        _adjacencyMap.put(vertex, new HashMap<>());
    }
    
    
    /**
     * Removes a vertex from a graph
     * 
     * @param vertex The vertex to be removed
     * @throws NoSuchElementException if the specified vertex
     * is not present in the graph
     */
    public void removeVertex(T vertex) throws NoSuchElementException 
    {
        // If the specified vertex is not present in the graph
        if (!_adjacencyMap.containsKey(vertex)) 
        {
            throw new NoSuchElementException("ERROR: Specified vertex"
                    + "is not present in graph.");
        }
        
        // For each vertex, if the vertex to be removed is a source or
        // destination vertex for any edge, remove the associated edge
        for (T otherVertex: _adjacencyMap.keySet())
        {
            // If the vertex to be removed is the destination node
            // for an edge, delete the edge
            if (edgeExists(otherVertex, vertex))
            {
                removeEdge(otherVertex, vertex);
            }
            // If the vertex to be removed is the source node for an edge,
            // delete the edge
            else if (edgeExists(vertex, otherVertex))
            {
                removeEdge(vertex, otherVertex);
            }
        }
        
        // Removes the vertex
        _adjacencyMap.remove(vertex);
    }
    
    
    /**
     * Adds an edge to the graph
     * 
     * @param source The source vertex of the edge
     * @param destination The destination vertex of the edge
     * @param weight The weight of the edge
     * @throws IllegalArgumentException if source and destination are the same,
     * the weight is negative, or an edge between the vertices already exists
     * @throws NoSuchElementException if the source or destination do not exist
     * in the graph
     */
    public void addEdge(T source, T destination, int weight) 
            throws IllegalArgumentException, NoSuchElementException 
    {
        // If either the source or the destination vertex do not exist
        // in the graph
        if (!_adjacencyMap.containsKey(source) 
                || !_adjacencyMap.containsKey(destination)) 
        {
            throw new NoSuchElementException("ERROR: Not all of the specified"
                    + "vertices are present in the Graph.");
        }
            
        // If the source and destination node are the same or the
        // edge weight is negative
        if (source.equals(destination) || weight < 0 )
        {
            throw new IllegalArgumentException("ERROR: Loops and negative edge "
                    + "weights not allowed.");
        }
        
        // Creates the edge from source to destination
        _adjacencyMap.get(source).put(destination, weight);
        
        // If the graph is undirected, creates connection from destination to
        // source as well with the same weight
        if (!_isDirected) 
        {
            _adjacencyMap.get(destination).put(source, weight);
        }
    }
    
    
    /**
     * Removes an edge from the graph
     * 
     * @param source The source vertex of the edge
     * @param destination The destination vertex of the edge
     * @throws NoSuchElementException if the source or the destination either
     * do not exist or an edge between the two specified vertices does not exist
     */
    public void removeEdge(T source, T destination) 
            throws NoSuchElementException 
    {
        if (!edgeExists(source, destination)) 
        {
            throw new NoSuchElementException("ERROR: No such edge exists.");
        }
        
        // Removes the specified edge
        _adjacencyMap.get(source).remove(destination);
        
        // If the graph is undirected or , removes the edge representing
        // the other direction
        if (!_isDirected)
        {
            _adjacencyMap.get(destination).remove(source);
        }
    }
    
    
    /**
     * Given a list representing a path, calculates and returns the length
     * of the path
     * 
     * @param path The sequence of vertices representing a path
     * @return The length of the path
     */
    public long pathLength(List<T> path)
    {   
        // The default value for the path length
        long pathLengthSum = UNDEFINED;
        boolean connected = true;
        int vertexIndex = 0;
        int pathSize = path.size();
        
        // Puts the path list into an ArrayList to ensure efficient operation
        ArrayList<T> pathList = new ArrayList<>(path);
        
        T sourceVertex;
        T destinationVertex;
        
        // If the path is not empty
        if (!pathList.isEmpty())
        {
            pathLengthSum = 0;
            
            // While the graph has not been verified to be unconnected and
            // the path list has not yet been exhausted
            while (connected && (vertexIndex < (pathSize - 1)))
            {
                // Sets the source and destination vertices for the edge of
                // interest during this iteration
                //
                // Note that the vertex index gets incremented after it is
                // referenced in the estinationVertex assignment
                sourceVertex = pathList.get(vertexIndex);
                destinationVertex = pathList.get(1 + vertexIndex++);
                
                // If the sequence is not a path
                if (_adjacencyMap.get(sourceVertex) == null ||
                        _adjacencyMap.get(sourceVertex).get(destinationVertex) 
                        == null)
                {
                    pathLengthSum = UNDEFINED;
                    connected = false;
                }
                else
                {
                    // Adds the edge weight to the running total
                    pathLengthSum += 
                            _adjacencyMap.get(sourceVertex)
                                    .get(destinationVertex);
                }
            }
        }
        
        return pathLengthSum;
    }
 
    
    /**
     * Based on a given source and destination vertex, can determine the
     * shortest path between the two
     * 
     * @param source The source vertex on the path
     * @param destination The destination vertex on the path
     * @return A list of the edges on the shortest path
     * @throws NoSuchElementException If one or both of the specified
     * vertices are not present in the graph
     */
    public List<Edge<T>> shortestPathBetween(T source, T destination)
            throws NoSuchElementException 
    {
        // If either the source or the destination are not present in the graph
        if (!_adjacencyMap.containsKey(source) 
                || !_adjacencyMap.containsKey(destination))
        {
            throw new NoSuchElementException("ERROR: One or both of the "
                    + "specified vertices does not exist in the graph.");
        }
        
        ArrayList<Edge<T>> shortestPath = new ArrayList<>();
        HashSet<T> unknownVertices = new HashSet<>(_adjacencyMap.keySet());
        PriorityQueue<Edge<T>> frontier = new PriorityQueue<>();
        HashMap<T, Integer> distance = new HashMap<>(_adjacencyMap.get(source));
        HashMap<T, T> predecessors = new HashMap<>();
        ArrayDeque<Edge<T>> edgeStack = new ArrayDeque<>();
        T vertexOfInterest;
        T lastVertex = destination;
        Edge<T> connectionOfInterest;
        int proposedLength;
        boolean firstIteration = true;
        
        // Populates the distance map with infinite distances and sets each
        // element's predecessor to null
        for (T vertex: _adjacencyMap.keySet())
        {
            distance.put(vertex, UNDEFINED);
            predecessors.put(vertex, null);
        }
        
        // While the frontier priority queue is not empty or this is
        // the first iteration through the while loop
        while (!frontier.isEmpty() || firstIteration)
        {   
            // If this is the first iteration through the while loop
            if (firstIteration)
            {
                // Set the distance of the source (0)
                distance.put(source, 0);
                
                // The source becomes the first vertexOfInterest
                vertexOfInterest = source;
                firstIteration = false;
            }
            else
            {
                // Selects the next connection of interest from the frontier
                connectionOfInterest = frontier.poll();
                vertexOfInterest = connectionOfInterest.getDestination();
            }

            // removes the vertex of interest from the unknown vertices set
            unknownVertices.remove(vertexOfInterest);

            // For each vertex adjacent to the vertex of interest
            for (T adjVert: _adjacencyMap.get(vertexOfInterest).keySet())
            {
                // If the adjacent vertex in question is currently unknown
                if (unknownVertices.contains(adjVert))
                {
                    // proposed length is the distance to the vertex of interest
                    // plus the edge length from the vertex of interest to
                    // the adjacent vertex
                    proposedLength = distance.get(vertexOfInterest)
                            + _adjacencyMap.get(vertexOfInterest).get(adjVert);
                    
                    // If the proposedLength is less than the current length
                    // of the path to adjVert
                    if ((proposedLength < distance.get(adjVert))
                            || (distance.get(adjVert) == UNDEFINED))
                    {
                        // Update the distance with the newly discovered
                        // shortest path
                        distance.put(adjVert, proposedLength);
                        predecessors.put(adjVert, vertexOfInterest);
                        
                        // Adds a new edge to the frontier such that
                        // the weight of the edge is the distance from
                        // the source to the vertex adjacent to the 
                        // vertex of interest
                        frontier.add(new Edge(vertexOfInterest, adjVert, 
                                proposedLength));
                    }
                }
            }
        }
        
        // If the specified source and destination vertices are not connected
        if (distance.get(destination) == UNDEFINED)
        {
            shortestPath = null;
        }
        // If the distance has been found
        else
        {
            // While the shortest path hasn't been completely loaded 
            // into the stack
            while (predecessors.get(lastVertex) != null)
            {
                // Pushes the edge between the currentVertex and its predecessor
                // into the stack
                edgeStack.push(new Edge<>(predecessors.get(lastVertex), 
                        lastVertex, _adjacencyMap.get(predecessors.
                                get(lastVertex)).get(lastVertex)));

                lastVertex = predecessors.get(lastVertex);
            }

            // pops all of the stack elements into the list in the desired order
            // yielding a shortest path from source to destination
            while (!edgeStack.isEmpty())
            {
                shortestPath.add(edgeStack.pop());
            }
        }
        
        return shortestPath;
    }
    
    
    /**
     * A private helper method that returns an arbitrary vertex in the graph
     * 
     * @return An arbitrary vertex from the graph
     */
    private T getArbitraryVertex()
    {
        T arbitraryVertex = null;
        int iteratorCount = 0;
        
        // Grabs an arbitrary vertex from the set
        for (T vertex: _adjacencyMap.keySet())
        {
            // If the iteratorCount matches the arbitrary value
            // the current vertex is the arbitrary vertex
            if (iteratorCount == 0)
            {
                arbitraryVertex = vertex;
            }
            
            iteratorCount++;
        }
        
        return arbitraryVertex;
    }
    
    
    /**
     * Returns a graph representing the minimum spanning tree
     * 
     * @return A graph representing the minimum spanning tree
     * @throws IllegalStateException If the graph is directed
     */
    public Graph<T> minimumSpanningTree() throws IllegalStateException 
    {   
        // If the graph is directed throw an IllegalStateException
        if (_isDirected)
        {
            throw new IllegalStateException("ERROR: Can only calculate "
                    + "minimum spanning tree of undirected graphs.");
        }
        
        // A flag indicating whether or not no span tree exists
        boolean noSpanTreeExists = false;
        boolean firstIteration = true;
        boolean validEdgeFound = false;
        int numVertices = _adjacencyMap.keySet().size();
        
        T vertexOfInterest = null;
        Edge<T> edgeOfInterest = null;
        
        // The minimum spanning tree
        Graph<T> minSpanTree = new Graph<>(_isDirected);
        
        PriorityQueue<Edge<T>> frontierHeap = new PriorityQueue<>();
        HashSet<T> unknownVertices = new HashSet<>(_adjacencyMap.keySet());
        
        // If the graph contains no vertices
        if (unknownVertices.isEmpty())
        {
            noSpanTreeExists = true;
            minSpanTree = null;
        }

        // For each vertex in the graph, check if it's an island,
        // if it is, then no span tree exists
        for (T vertex: _adjacencyMap.keySet())
        {
            // if the graph is not a singleton
            if (numVertices != 1)
            {
                // If the vertex is an island and not a singleton
                if (_adjacencyMap.get(vertex).keySet().isEmpty())
                {
                    noSpanTreeExists = true;
                    //
                    minSpanTree = null;
                }
            }
            else
            {
                // If the tree is a single vertex
                minSpanTree.addVertex(vertex);
                unknownVertices.remove(vertex);
            }
        }
        
        // By use of PRIM'S ALGORITHM
        //
        // While there still exists unknown vertices
        while(!unknownVertices.isEmpty() && !noSpanTreeExists)
        {
            // If this is the first iteration through the while loop
            if (firstIteration)
            {
                // Initializes the first vertex of interest to an
                // arbitrary vertex within the graph
                vertexOfInterest = getArbitraryVertex();
                
                // Adds the vertex to the minimum spanning tree
                minSpanTree.addVertex(vertexOfInterest);
                firstIteration = false;
            }
            else
            {
                // While an edge from the connected graph that connects to
                // an unknown vertex is not found and it has not been verified
                // that no span tree exists
                while (!validEdgeFound && !noSpanTreeExists)
                {
                    // The next edge of interest is the edge in the frontier 
                    // with the smallest weight value
                    edgeOfInterest = frontierHeap.poll();
                    
                    // If the destination vertex is unknown and the frontier 
                    // heap is not empty or if the frontier is empty and 
                    // there is only one more unknown vertex
                    if (!frontierHeap.isEmpty() 
                            && unknownVertices.
                                    contains(edgeOfInterest.getDestination())
                            || (frontierHeap.isEmpty() 
                            && unknownVertices.size() == 1))
                    {
                        validEdgeFound = true;
                    }
                    // If the frontier heap is empty and there are more
                    // than 1 unknown vertices then we know that we have
                    // an unconnected graph
                    else if (frontierHeap.isEmpty() 
                            && unknownVertices.size() > 1)
                    {
                        noSpanTreeExists = true;
                        minSpanTree = null;
                    }
                }
                
                if ((vertexOfInterest = edgeOfInterest.getDestination()) 
                        != null && !noSpanTreeExists)
                {
                    // Adds the vertex of interest and the new edge to 
                    // the minimum spanning tree
                    minSpanTree.addVertex(vertexOfInterest);
                    minSpanTree.addEdge(edgeOfInterest.getSource(), 
                            vertexOfInterest, edgeOfInterest.getWeight());
                }
            }
            
            if (!noSpanTreeExists)
            {
                // Removes the new vertex of interest from the frontier heap
                unknownVertices.remove(vertexOfInterest);

                // For each vertex adjacent to the vertexOfInterest
                for (T vertex: _adjacencyMap.get(vertexOfInterest).keySet())
                {
                    // If the adjacent vertex is not yet in the minimum spanning
                    // tree, put the edge between it
                    // and the vertex of interest in the heap
                    if (unknownVertices.contains(vertex)) { 
                        frontierHeap.add(new Edge<>(vertexOfInterest, vertex, 
                                _adjacencyMap.get(vertexOfInterest)
                                        .get(vertex)));
                    }
                }

                // Resets the validEdgeFound flag
                validEdgeFound = false;
            }
        }
        
        return minSpanTree;
    }
    
    
    /**
     * Returns a tour as a list of vertices that is very close to optimal
     * using the Inver-over operator and pseudocode from "Inver-over Operator 
     * for the TSP" by Guo Tao and Zbigniew Michalewicz that can be found 
     * at the URL below.
     * 
     *      http://dl.acm.org/citation.cfm?id=668606
     * 
     * @param populationSize The size of the population
     * @param inversionProbability The probability of generating a 
     * random inversion
     * @param terminationIterations The number of iterations before termination
     * @return A list of vertices representing a close-to-optimal tour 
     */
    public List<T> getOptimalTour(int populationSize, 
            float inversionProbability, int terminationIterations)
    {
        ArrayList<T> vertexList = new ArrayList<>(getVertices());
        ArrayList<T> permutation = new ArrayList<>();
        ArrayList<T> currentTour = new ArrayList<>();
        ArrayList<T> otherTour = new ArrayList<>();
        ArrayList<T> tour = new ArrayList<>();
        ArrayList<T> previousSmallest = new ArrayList<>();
        ArrayList<ArrayList<T>> population = new ArrayList<>();
        ArrayList<T> smallest = new ArrayList<>();
        Random randomGenerator = new Random();
        T vertex;
        T otherVertex;
        int vertexRandomIndex;
        int otherVertexIndex;
        int offset = 0;
        long currentTourLength = 0;
        long tourLength = 0;
        long previousSmallestLength = 0;
        long smallestLength = 0;
        int iterationCount = terminationIterations;
        int tourSize = vertexList.size();
        boolean vertexAdjacency = false;
        
        // Add a permutation of the vertex set to the population 
        // populationSize times
        for (int i = 0; i < populationSize; i++)
        {
            // Creates a copy of the vertex list that will be shuffled
            // and added to the population
            permutation = (ArrayList<T>) vertexList.clone();
            
            // Shuffles vertexSet into a permutation of itself
            Collections.shuffle(permutation);
            
            // Adds this permutation to the population
            population.add(permutation);
        }
        
        // While the iteration count is not zero, keep looking for
        // a new optimal tour
        while (iterationCount != 0)
        {
            // Initializes the best tour in population, will be updated
            // as iteration through the population commences
            smallest = (ArrayList<T>) population.get(0).clone();
            
            // Initializes the length of the aforementioned initial best 
            // tour in population
            smallestLength = pathLength(smallest) 
                    + getEdgeWeight(smallest
                            .get(tourSize - 1), smallest.get(0));
            
            // For each member of the population
            for (int i = 0; i < populationSize; i++)
            {
                // Tour is set as the ith member of the population
                tour = population.get(i);
                
                // Calculates the length of the tour including the final
                // edge length
                tourLength = pathLength(tour) 
                        + getEdgeWeight(tour.get(tourSize - 1), tour.get(0));

                // Sets the current tour on which to operate as the current
                // smallest tour discovered
                currentTour = (ArrayList<T>) tour.clone();

                // Randomly selects an index corresponding to a vertex in 
                // the current tour
                vertexRandomIndex = randomGenerator.nextInt(tourSize);

                // Randomly selects a currentVertex from currentTour
                vertex = currentTour.get(vertexRandomIndex);

                // While vertex and otherVertex are not adjacent
                // perform inver-over operation on the current tour
                while (!vertexAdjacency)
                {
                    // If the inversionProbability is greater than or equal
                    // to the a randomly generated float value [0...1], 
                    // select an otherVertex from the remaining cities
                    // in currentTour
                    if (randomGenerator.nextFloat() <= inversionProbability)
                    {
                        // Set otherVertex as an element in currentTour
                        // other than vertex
                        do
                        {
                            otherVertex = currentTour
                                    .get(randomGenerator.nextInt(tourSize));
                        }
                        while (otherVertex.equals(vertex));
                    }
                    // If the inversion probability is less than the randomly
                    // generated float value [0...1], then select 
                    // another tour from the population and assign otherVertex
                    // to the element following vertex in the selected tour
                    else
                    {
                        // Assign otherTour to a member of the population
                        // other than the currentTour
                        do
                        {
                            otherTour = population.get(randomGenerator
                                    .nextInt(populationSize));
                        }
                        while (otherTour.equals(currentTour));

                        // Assign the element following vertex in the other
                        // tour to otherVertex; It is worth noting that if 
                        // vertex is the last element in otherTour, otherVertex
                        // is assigned to the first element in otherTour since a 
                        // tour is cyclic
                        otherVertex = otherTour.get((otherTour
                                .indexOf(vertex) + 1) % tourSize);
                    }

                    // If the element before or after vertex in the currentTour
                    // is equal to otherVertex, then there is vertexAdjacency
                    if (currentTour
                            .get((vertexRandomIndex + 1) % tourSize)
                            .equals(otherVertex)
                            || currentTour
                                    .get(Math.abs((vertexRandomIndex - 1) 
                                            % tourSize))
                                    .equals(otherVertex))
                    {
                        vertexAdjacency = true;
                    }

                    else
                    {
                        // The index in currentTour that holds otherVertex
                        otherVertexIndex = currentTour.indexOf(otherVertex);

                        // If the index of otherVertex is less than the index
                        // of vertex, reverse the sublist from otherVertex
                        // to vertex in the currentTour
                        if (otherVertexIndex < vertexRandomIndex)
                        {
                            // Reverses the order of the sublist from
                            // otherVertex to vertex in currentTour
                            offset = (tourSize - vertexRandomIndex);
                            
                            // Rotate the current tour so that vertex becomes
                            // the zeroth index of the currentTour
                            Collections.rotate(currentTour, offset);
                            
                            // In the rotated current tour, reverse the sublist
                            // from the element after vertex (at the zeroth
                            // index) to the new position of otherVertex
                            Collections.reverse(currentTour
                                    .subList(1, 
                                            (otherVertexIndex + offset) + 1));
                            
                            // Rotates the tour back into its initial relative
                            // configuration
                            Collections.rotate(currentTour, vertexRandomIndex);
                            
                            otherVertexIndex = (vertexRandomIndex + 1) 
                                    % tourSize;
                        }
                        // If the index of vertex is less than the index of
                        // otherVertex, reverse the sublist from vertex
                        // to otherVertex in the currentTour
                        else
                        {
                            // Reverses the order of the sublist from the 
                            // element after vertex to otherVertex in 
                            // currentTour
                            Collections.reverse(currentTour
                                    .subList((vertexRandomIndex + 1), 
                                            otherVertexIndex + 1)); 
                            otherVertexIndex = (vertexRandomIndex + 1) 
                                    % tourSize;
                        }

                        vertex = otherVertex;
                        vertexRandomIndex = otherVertexIndex;
                    }
                }

                // Resets the while loop flag
                vertexAdjacency = false;

                // Calculates the weight of the final edge in the current
                // tour 
                currentTourLength = pathLength(currentTour) 
                        + getEdgeWeight(currentTour.get(tourSize - 1), 
                                currentTour.get(0));
                
                // If the path length of the currentTour is smaller than or
                // equal to the path length of the smallest path, the
                // currentTour becomes the new tour
                if (currentTourLength <= tourLength)
                {
                    population.set(i, currentTour);
                    
                    // If the currentTour has a smaller path length than the
                    // smallest path length, make it the new smallest
                    if (currentTourLength < smallestLength)
                    {
                        smallest = currentTour;
                        smallestLength = currentTourLength;
                    }
                }
                // If tour has a smaller path length than the smallest path 
                // length, make it the new smalelst
                else if (tourLength < smallestLength)
                {
                    smallest = tour;
                    smallestLength = tourLength;
                }
            }
            
            // If the shortestPathLength remains unchanged decrement the 
            // iterationCount; If previousSmallest is empty, then this is 
            // the first loop iteration in which case, iterationCount may
            // be decremented
            if (smallestLength == previousSmallestLength 
                    || previousSmallest.isEmpty())
            {
                iterationCount--;
            }
            else 
            {
                iterationCount = terminationIterations;
            }
            
            // Update the previousSmallest to the the tour from
            // this iteration
            previousSmallest = smallest;
            previousSmallestLength = smallestLength;
        }
        
        return smallest;
    }
    
    
    /**
     * A factory method that reads graph information from a .csv file and
     * generates an instance of the specified graph
     * 
     *      Reads in values from .csv files assuming the following format:
     *          <# of vertices>
     *          <first vertex>
     *          <second vertex>
     *          ...
     *          <nth vertex>
     *          <# of edges>
     *          <source>,<destination>,<weight>
     *          <source>,<destination>,<weight>
     *          ...
     * 
     * @param isDirected A boolean value that is true if the graph is directed
     * @param inputFile The .csv file being read from
     * @return The graph described in the .csv file
     * @throws IOException If the file format is incorrect
     */
    public static Graph<String> fromCSVFile(boolean isDirected, 
            Scanner inputFile) throws IOException
    {
        final String ERROR = "ERROR: Please ensure that the .csv file "
                + "containing the graph is formatted correctly.";
        
        // Will be the starting index when parsing edge weights
        final int START_INDEX = 1;
        
        // The graph being generated
        Graph<String> csvGraph = new Graph<>(isDirected);
        
        int numberOfVertices;
        int numberOfEdges;
        
        try {
            // Scrapes the number of vertices from the first line of the file
            numberOfVertices = Integer.parseInt(inputFile.nextLine().trim());
        }
        catch (NumberFormatException nfe)
        {
            throw new IOException(ERROR);
        }
        
        // Adds all of the specifed vertices to the graph
        for (int i = 0; i < numberOfVertices; i++)
        {
            try 
            {
                csvGraph.addVertex(inputFile.nextLine().trim());
            }
            catch (IllegalArgumentException iae)
            {
                throw new IOException(ERROR);
            }
        }
        
        try {
            // Scrapes the number of edges from the (N + 2)nd line of the file
            numberOfEdges = Integer.parseInt(inputFile.nextLine().trim());
        }
        catch (NumberFormatException nfe)
        {
            throw new IOException(ERROR);
        }
        
        // Changes the scanners delimiter to ","
        inputFile.useDelimiter(",");
        
        // Adds all of the specified edges to the graph
        for (int i = 0; i < numberOfEdges; i++)
        {
            try {
                // Parses edge information from edge specification lines in
                // the .csv file
                // Note: The substring is utilized because the edge weight
                // String returned by the scanner starts with a delimiting
                // comma left over from the last scanner action
                csvGraph.addEdge(inputFile.next().trim(), 
                        inputFile.next().trim(), 
                        Integer.parseInt(inputFile.nextLine()
                                .substring(START_INDEX).trim()));
            }
            catch(NumberFormatException nfe)
            {
                throw new IOException(ERROR);
            }
            catch(NoSuchElementException | IllegalArgumentException nsee)
            {
                throw new IOException(ERROR);
            }
        }
        
        // Returns the generated graph
        return csvGraph;
    }
    
    
    /**
    * Construct an undirected graph from an XML encoded TSP file
    *
    * @param inputFile an XML encoded TSP file from 
    *        <a href="http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/">
    *        http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/</a>
    *
    * @return graph populated from the file
    * 
    * @throws ParserConfigurationException, SAXException, IOException if
    *         the file doesn't conform to the specification
    */
    public static Graph<String> fromTSPFile(InputStream inputFile) 
         throws ParserConfigurationException, SAXException, IOException {

        /**
         * The Handler for SAX Parser Events.
         *
         * This inner-class extends the default handler for the SAX parser
         * to construct a graph from a TSP file
         *
         * @see org.xml.sax.helpers.DefaultHandler
         */
        class TSPGraphHandler extends DefaultHandler {
            // Instantiate an undirected graph to populate; vertices are
            // integers though we treat them as strings for extension to other
            // similarly-formed files representing, say, GFU.
            private Graph<String> _theGraph = new Graph<>(false);

            private final int NO_WEIGHT = -1;
            // As we parse we need to keep track of when we've seen
            // vertices and edges
            private int _sourceVertexNumber = 0;
            private String _destinationVertexName = null;
            private String _sourceVertexName = null;
            private int _edgeWeight = NO_WEIGHT;
            private boolean _inEdge = false;


            /**
             * Parser has seen an opening tag
             *
             * For a <pre>vertex</pre> tag we add the vertex to the graph
             * the first time we encounter it.
             * For an <pre>edge</pre> tag we remember the weight of the edge.
             *
             * {@inheritDoc}
             */
            @Override
            public void startElement(String uri, String localName,
                   String qName, Attributes attributes) throws SAXException {

                // We only care about vertex and edge elements
                switch (qName) {

                    case "vertex":
                        // See if the vertices are named; if so, use the
                        // name, otherwise use the number
                        _sourceVertexName = attributes.getValue("name");
                        if (_sourceVertexName == null) {
                            _sourceVertexName = Integer.toString(_sourceVertexNumber);
                        }
                        // If is vertex 0 then it's the first time we're seeing it; 
                        // add it to the graph. Other vertices will be added
                        // as we encounter their edges
                        if (_sourceVertexNumber == 0) {
                            _theGraph.addVertex(_sourceVertexName);
                        }
                        break;

                    case "edge":
                        // Edges have the destination vertex within so
                        // indicate that we're inside an edge so that the
                        // character-parsing method below will grab the
                        // destination vertex as it encounters it
                        _inEdge = true;
                        // The weight of the edge is given by the "cost"
                        // attribute
                        _edgeWeight = (int) Double.parseDouble(attributes.getValue("cost"));
                        break;

                    default: // ignore any other opening tag
                }
            }


            /**
             * Parser has seen a closing tag.
             *
             * For a <pre>vertex</pre> tag we increment the vertex number
             * to keep track of which vertex we're parsing.
             * For a <pre>edge</pre> tag we use the number of the edge and
             * the weight we saw in the opening tag to add an edge to the
             * graph.
             *
             * {@inheritDoc}
             */
            @Override
            public void endElement(String uri, String localName,
                   String qName) throws SAXException {

                // Again, we only care about vertex and edge tags
                switch (qName) {

                    case "vertex":
                        // End of a vertex so we're moving on to the next
                        // source vertex number
                        _sourceVertexNumber++;
                        // Clear out the name so we don't inherit it in some
                        // mal-formed entry later
                        _sourceVertexName = null;
                        break;

                    case "edge":
                        // We've finished an edge so we have collected all the
                        // information needed to add an edge to the graph
                        _inEdge = false;
                        // If this is the first set of edges (i.e., we're on
                        // the first source vertex) then this is the first
                        // time we've seen the destination vertex; add it to
                        // the graph
                        if (_sourceVertexNumber == 0) {
                            _theGraph.addVertex(_destinationVertexName);
                        }
                        // Should now be safe to add an edge between the
                        // source and destination
                        _theGraph.addEdge(_sourceVertexName, 
                                        _destinationVertexName, _edgeWeight);
                        // Clear out the attributes of this edge so we don't
                        // accidentally inherit them should we parse a
                        // mal-formed edge entry later
                        _destinationVertexName = null;
                        _edgeWeight = NO_WEIGHT;
                        break;

                    default: // ignore any other closing tag
                }
            }


            /**
            * Parser has seen a string of characters between opening and
            * closing tag. The only characters we care about occur within
            * an <pre>edge</pre> tag and are the destination vertex.
            *
            * {@inheritDoc}
            */
            @Override
            public void characters(char[] ch, int start, int length) throws SAXException {
                // If we're within an edge, then this string of characters
                // is the number of the destination vertex for this edge.
                // Remember the destination vertex
                if (_inEdge) {
                    _destinationVertexName = new String(ch, start, length);
                }
            }


            /**
            * @return the graph constructed
            */
            Graph<String> getGraph() {
                return _theGraph;
            }

        } // TSPHandler


        // Create a handler and use it for parsing
        TSPGraphHandler tspHandler = new TSPGraphHandler();

        // Here's where we do the actual parsing using the local class
        // defined above. Give the parser an instance of the class above 
        // as the handler and parse away!
        SAXParserFactory.newInstance().newSAXParser().parse(inputFile, tspHandler);

        // Graph should now be populated, return it
        return tspHandler.getGraph();

    }
    
    
    /**
     * Returns a string containing a human-readable representation of a graph
     * 
     * e.g. For a directed graph containing vertices A, B, and C with edges 
     * (A, B), (A, C), and (B, C), (all with weight 0) the output
     * would be as follows:
     * 
     *               DIRECTED [A={B=0, C=0}, B={C=0}, C={}]
     * 
     * for an undirected graph with the same vertices and edge,
     * 
     *        UNDIRECTED [A={B=0, C=0}, B={A=0, C=0}, C={A=0, B=0}]
     * 
     * @return A string containing a human-readable representation of a graph
     */
    @Override
    public String toString() 
    {
        String isDirected;
        
        // Sets isDirected as appropriate
        if (_isDirected)
        {
            isDirected = "DIRECTED ";
        }
        else
        {
            isDirected = "UNDIRECTED ";
        }
        
        // Returns a string that contains the graph's vertices and edges
        return isDirected + _adjacencyMap.toString();
    }
    
    
    /**
     * Represents an Edge in a Graph
     * 
     * @param <E> The type of the source and destination vertices
     */
    public static class Edge<E> implements Comparable<Edge<E>>
    {
        private E _source;
        private E _destination;
        private int _weight;
        
        
        /**
         * Constructs an edge with a specified source, destination, and weight
         * 
         * @param source The source vertex
         * @param destination The destination vertex
         * @param weight The weight of the edge
         */
        public Edge(E source, E destination, int weight) 
        {
            _source = source;
            _destination = destination;
            _weight = weight;
        }
        
        
        /**
         * Returns the weight of the edge
         * 
         * @return The edge's weight
         */
        public int getWeight() 
        {
            return _weight;
        }
        
        
        /**
         * Returns the edge's source vertex
         * 
         * @return The edge's source vertex
         */
        public E getSource() 
        {
            return _source;
        }
        
        
        /**
         * Returns the edge's destination vertex
         * 
         * @return The edge's destination vertex
         */
        public E getDestination()
        {
            return _destination;
        }

        
        /**
         * Compares one edge to another
         * 
         * @param e The edge to be compared with this edge instance
         * @return A positive integer if this instance possesses a greater 
         * edge weight than e, 0 if the edge weights are equal, or a negative
         * integer if this instance's edge weight is less than the edge
         * weight of e
         */
        @Override
        public int compareTo(Edge<E> e)
        {
            return _weight - e.getWeight();
        }

        
        /**
         * Returns a human-readable representation of the edge as a String
         * 
         * e.g. For an edge with source vertex, "hello", and destination vertex,
         * "world", with a weight of 10, this method would return:
         * 
         *          hello___10___world
         * 
         * @return A human-readable String representation of the edge
         */
        @Override
        public String toString()
        {
            return _source.toString() + "___" + _weight + "___" 
                    + _destination.toString();
        }
    }
}