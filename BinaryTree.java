import java.util.NoSuchElementException;

/**
 * Represents a node in a binary tree.
 * 
 * @author Thomas Noel
 * @param <T> The data type being represented in the tree
 */
public class BinaryTree <T extends Comparable<? super T>> 
{
    
    private BinaryNode<T> _root;
    
    /**
     * Constructs an empty binary tree
     */
    public BinaryTree()
    {
        _root = null;
    }
    
    
    /**
     * Empties the binary tree of all data
     */
    public void makeEmpty()
    {
        _root = null;
    }
    
    
    /**
     * Checks to see if tree is empty
     * 
     * @return True if tree is empty, false otherwise
     */
    public boolean isEmpty()
    {
        return _root == null;
    }
    
    
    /**
     * Checks to see if element x is contained in the tree
     * 
     * @param x The values whose presence is being determined
     * @return True if element present, else false
     */
    public boolean contains(T x)
    {
        // TODO: Implement contains
        return contains(x, _root);
    }
    
    
    /**
     * Finds the element with minimum value in the tree.
     * 
     * @return The minimum element in the tree
     */
    public T findMin()
    {
        if (isEmpty())
        {
            throw new NoSuchElementException();
        }
       
        return findMin(_root).element;
    }
    
    
    /**
     * Finds the element with maximum value in the tree.
     * 
     * @return The maximum element in the tree
     */
    public T findMax()
    {
        if (isEmpty())
        {
            throw new NoSuchElementException();
        }

        return findMax(_root).element;
    }
    
    
    /**
     * Inserts element x into the tree
     * 
     * @param x The data to be assimilated into the tree structure
     */
    public void insert(T x)
    {
        _root = insert(x, _root);
    }
    
    
    /**
     * Removes the specified element from the
     * binary tree
     * 
     * @param x The element being removed
     */
    public void remove(T x)
    {
        _root = remove(x, _root);
    }

    
    /**
     * Prints the tree structure intelligibly
     */
    public void printTree()
    {
        printTree(_root);
    }
    
    
    /**
     * Internal method to find an item in a subtree
     * 
     * @param x The item sought
     * @param t The root of the subtree
     * @return True if item is present, false otherwise
     */
    private boolean contains(T x, BinaryNode<T> t)
    {
        int comparisonResult;
        boolean containsElement = false;
        
        if (x != null && t != null)
        {
            comparisonResult = t.element.compareTo(x);
            
            // If the element in binary node t is equal to x
            if (comparisonResult == 0)
            {
                containsElement = true;
            }
            else if (comparisonResult > 0)
            {
                containsElement = contains(x, t.left);
            }
            else if (comparisonResult < 0)
            {
                containsElement = contains(x, t.right);
            }
        }
        
        return containsElement;
    }
    
    
    /**
     * Finds the minimum element of a binary search tree
     * 
     * @param t The subtree root
     * @return The node of the subtree containing the min-
     * value element
     */
    private BinaryNode<T> findMin(BinaryNode<T> t)
    {
        BinaryNode<T> minNode = t;

        if (t.left != null)
        {
            minNode = findMin(t.left);
        }
        
        return minNode;
    }
    
    
    /**
     * Finds the maximum element of a binary search tree
     * 
     * @param t The subtree root
     * @return The node of the subtree containing the max-
     * value element
     */
    private BinaryNode<T> findMax(BinaryNode<T> t)
    {
        BinaryNode<T> maxNode = t;

        if (t.right != null)
        {
            maxNode = findMax(t.right);
        }
        
        return maxNode;
    }
    
    
    /**
     * Inserts the specified element into the binary search tree
     * using appropriate binary search tree structure semantics
     * 
     * @param x The element to be inserted
     * @param t The root of the subtree
     * @return The node containing the inserted element
     */
    private BinaryNode<T> insert(T x, BinaryNode<T> t)
    {
        BinaryNode<T> rootNode = t;
        int comparisonValue;
        
        if (x != null)
        {   
            if (t == null)
            {
                rootNode = new BinaryNode(x);
            }
            else
            {
                comparisonValue = t.element.compareTo(x);

                // If the value of the element to be inserted
                // is equal to the current subtree root's element
                if (comparisonValue == 0)
                {
                    t.numInstances++;
                }
                // If the value of x is greater than the value
                // of the element in node t
                else if (comparisonValue < 0)
                {
                    // go right in tree
                    if (t.right == null)
                    {
                        t.right = new BinaryNode(x);
                    }
                    else
                    {
                        insert(x, t.right);
                    }
                }
                // If the value of x is less than the value
                // of the element in node t
                else if (comparisonValue > 0)
                {
                    // go left in tree
                    if (t.left == null)
                    {
                        t.left = new BinaryNode(x);
                    }
                    else
                    {
                        insert(x, t.left);
                    }
                }
            }
        }
        
        return rootNode;
    }
    
    
    /**
     * Removes the specified element from the binary tree if the element
     * is present.
     * 
     * @param x The element being inserted into the tree
     * @param t The root of the current subtree
     * @return The root of the modified tree
     */
    private BinaryNode<T> remove(T x, BinaryNode<T> t)
    {
        BinaryNode<T> rootNode = t;
        BinaryNode<T> newRoot;
        int comparisonResult;
        
        if (x != null && t != null)
        {
            comparisonResult = t.element.compareTo(x);
            
            // If the element to be removed has been
            // found, grab its replacement, sever its ties
            // with its parent node, and set its right and left
            // children to the right and left children of the current
            // root
            if (comparisonResult == 0)
            {
                
                if (t.numInstances == 1)
                {
                    // Ideally, the current node is replaced
                    // with the minimum child from the right subtree
                    if (t.right != null)
                    {
                        newRoot = removeMin(t.right);
                        newRoot.left = t.left;
                    }
                    // If no right subtree exists, then this node can be
                    // replaced with the maximum child from the left subtree
                    else if (t.left != null)
                    {
                        newRoot = removeMax(t.left);
                        newRoot.right = t.right;
                    }
                    // If this is a leaf
                    else
                    {
                        newRoot = null;
                    }

                    rootNode = newRoot;
                }
                
                t.numInstances--;
            }
            else if (comparisonResult < 0)
            {
                newRoot = remove(x, t.right);
                t.right = newRoot;
                
            }
            else if (comparisonResult > 0)
            {
                newRoot = remove(x, t.left);
                t.left = newRoot;
            }
        }
        
        return rootNode;
    }
    
    
    /**
     * Removes the smallest element of the subtree and severs it 
     * from its parent; A private helper method for remove().
     * 
     * @param t The root of the subtree
     * @return The minimum element of the modified tree
     */
    private BinaryNode<T> removeMin(BinaryNode<T> t)
    {
        BinaryNode<T> newNode;
        
        if (t.left != null)
        {
            newNode = removeMin(t.left);
            
            // Break the link if two levels away from minimum
            if (t.left.left == null)
            {
                // If the leftmost node has a right child, make this node
                // the left child of the root node (Special case)
                if (t.left.right != null)
                {
                    t.left = t.left.right;
                }
                else
                {
                    t.left = null;
                }
            }
        }
        else
        {
            newNode = t;
        }
        
        return newNode;
    }
    
    
    /**
     * Removes the largest element of the subtree and severs it 
     * from nodes attached; A private helper method for remove().
     * 
     * @param t The root of the subtree
     * @return The root of the modified tree
     */
    private BinaryNode<T> removeMax(BinaryNode<T> t)
    {
        BinaryNode<T> newNode;
        
        if (t.right != null)
        {
            newNode = removeMax(t.right);
            
            // Break the link if two levels away from minimum
            if (t.right.right == null)
            {
                // If the leftmost node has a right child, make this node
                // the left child of the root node (Special case)
                if (t.right.left != null)
                {
                    t.right = t.right.left;
                }
                else
                {
                    t.right = null;
                }
            }
        }
        else
        {
            newNode = t;
        }
        
        return newNode;
    }
    
    /**
     * Prints out the contents of the binary tree
     * 
     * @param t The root of the current subtree
     */
    private void printTree(BinaryNode<T> t)
    {
        // Depth first traversal
        System.out.print("Element: ");
        System.out.println(t.element);
        
        if (t.left != null)
        {
            printTree(t.left);
        }
        if (t.right != null)
        {
            printTree(t.right);
        }
    }
    
    
    /**
     * The nodes of the Binary Tree. Accommodates the
     * use of generic types.
     * 
     * @param <E> The data type of the node
     */
    private static class BinaryNode <E> 
    {
        
        int numInstances;
        E element;
        BinaryNode<E> left;
        BinaryNode<E> right;
        
        /**
         * Builds a new Binary Node
         * 
         * @param newElement The element to be stored in the node
         */
        public BinaryNode(E newElement)
        {
            this(newElement, null, null);
        }
        
        /**
         * Builds a new Binary node with a specified left and
         * specified right sibling
         * 
         * @param newElement The node's element
         * @param lt The left sibling of the node
         * @param rt The right sibling of the node
         */
        public BinaryNode(E newElement, BinaryNode<E> lt, BinaryNode<E> rt)
        {
            element = newElement;
            left = lt;
            right = rt;
            numInstances = 1;
        }
    }
    
}
