/**
// code is far away from bug with Buddha protection
//
//
//                       _oo0oo_
//                      o8888888o
//                      88" . "88
//                      (| -_- |)
//                      0\  =  /0
//                    ___/`---'\___
//                  .' \\|     |// '.
//                 / \\|||  :  |||// \
//                / _||||| -:- |||||- \
//               |   | \\\  -  /// |   |
//               | \_|  ''\---/''  |_/ |
//               \  .-\__  '-'  ___/-. /
//             ___'. .'  /--.--\  `. .'___
//          ."" '<  `.___\_<|>_/___.' >' "".
//         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
//         \  \ `_.   \_ __\ /__ _/   .-` /  /
//     =====`-.____`.___ \_____/___.-`___.-'=====
//                       `=---='
//
//
//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @author Wang YanHao
 */

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.PriorityQueue;
import java.util.Stack;

// write your matric number here: A0113742
// write your name here: Wang YanHao
// write list of collaborators here (reading someone's post in Facebook group and using the idea is counted as collaborating): null

class Supermarket {
	private int N; // number of items in the supermarket. V = N+1 
	private int K; // the number of items that Steven has to buy
	private int[] shoppingList; // indices of items that Steven has to buy
	private int[][] T; // the complete weighted graph that measures the direct walking time to go from one point to another point in seconds

	private int[][] memo;
	private PriorityQueue<IntegerPair> pq = new PriorityQueue<IntegerPair>();
	private Stack<Node> stackBAndB = new Stack<Node>();
	public static int[][] reducedT;
	private int size;
	private int[] reducedV;
	private int[] dist;
	private int mask;
	private int bestTour;

	public Supermarket() {}

	int Query() {
		if (K <= 18) {
			// dp
			memo = new int[K+1][(int) Math.pow(2, K+1)];
			mask = (int) (Math.pow(2, size) - 1);

			int ans = DP(0, 0);
			return ans;
		} else {
			bestTour = Integer.MAX_VALUE;

			// BranchAndBound
			Node root = new Node(size);

			for (int i=0; i<size; i++) {
				root.assignRowIndex(i, (byte)i);
				root.assignColIndex(i, (byte)i);

				// update reducedT entries
				reducedT[i][i] = -1;
			}

			root.assignMatrixValue(reducedT);
			root.reduceRow();
			root.reduceCol();

			stackBAndB.push(root);
			while (!stackBAndB.isEmpty()) {
				Node processedNode = stackBAndB.pop();

				if (processedNode.getSize() == 2) {
					// have a feasible solution
					processedNode.calFeasible();
					
					// update bestTour
					if (processedNode.getLowerBound() < bestTour) {
						bestTour = processedNode.getLowerBound();
					}

				} else if (processedNode.getLowerBound() < bestTour) {
					processedNode.selectNode();
					
					// branch
					
					// right branch : include entry with highest penalty
					Node rightChild = new Node(processedNode, processedNode.getHighestPenaltyRowIndex(), processedNode.getHighestPenaltyColIndex());
					// update LB
					rightChild.reduceRow();
					rightChild.reduceCol();

					// left branch : exclude entry with highest penalty
					// manipulate processedNode to get leftChild
					boolean isContinued = processedNode.alterForLeftChild();
					if (isContinued) {
						Node leftChild = processedNode;
					
						// check for LB for left and right child
						if (rightChild.getLowerBound() <= leftChild.getLowerBound()) {
							stackBAndB.push(leftChild);
							stackBAndB.push(rightChild);
						} else {
							stackBAndB.push(rightChild);
							stackBAndB.push(leftChild);
						}
					} else {
						stackBAndB.push(rightChild);
					}
				}
			}

			return bestTour;
		}
	}

	private int DP(int i, int j) {
		if (j == mask) {
			return reducedT[i][0];
		}
		if (memo[i][j] != 0) {
			return memo[i][j];
		}

		memo[i][j] = Integer.MAX_VALUE;
		for (int k=0; k<=K; k++) {
			int masked = j | (1<<k);
			if (masked != j ) {
				memo[i][j] = Math.min(memo[i][j], 
						reducedT[i][k] + DP(k, masked));
			}
		}
		return memo[i][j];
	}

	/**
	 * Construct reducedT
	 */
	private void preprocessing() {
		// initialize
		size = K+1;
		reducedT = new int[size][size]; 
		reducedV = new int[size];
		dist = new int[N+1];

		// construct reducedV
		reducedV[0] = 0;
		for (int i=0; i<K; i++) {
			reducedV[i+1] = shoppingList[i];
		}

		// construct reducedT
		for (int i=0; i<K; i++) {
			dijkstra(reducedV[i]);

			for (int j=i+1; j<=K; j++) {
				reducedT[i][j] = reducedT[j][i] = dist[reducedV[j]];
			}
		}
	}

	/**
	 * Using dijkstra to find shortest path from source vertex to all other 
	 * required to visited vertices in the shoppingList.
	 */
	private void dijkstra(int sv) {
		Arrays.fill(dist, Integer.MAX_VALUE);
		pq.add(new IntegerPair(sv, 0));
		dist[sv] = 0;

		while (!pq.isEmpty()) {
			IntegerPair temp = pq.poll();
			int tempDist = temp.second();
			int tempVertex = temp.first();

			if (tempDist <= dist[tempVertex]) {
				for (int i=0; i<=N; i++) {
					if (dist[i] > T[tempVertex][i] + dist[tempVertex]) {
						dist[i] = T[tempVertex][i] + dist[tempVertex];
						pq.add(new IntegerPair(i, dist[i]));
					}
				}
			}
		}
	}

	void run() throws Exception {
		// do not alter this method to standardize the I/O speed (this is already very fast)
		IntegerScanner sc = new IntegerScanner(System.in);
		PrintWriter pw = new PrintWriter(new BufferedWriter(new OutputStreamWriter(System.out)));

		int TC = sc.nextInt(); // there will be several test cases
		while (TC-- > 0) {
			// read the information of the complete graph with N+1 vertices
			N = sc.nextInt(); K = sc.nextInt(); // K is the number of items to be bought

			shoppingList = new int[K];

			for (int i = 0; i < K; i++)
				shoppingList[i] = sc.nextInt();

			T = new int[N+1][N+1];
			for (int i = 0; i <= N; i++)
				for (int j = 0; j <= N; j++)
					T[i][j] = sc.nextInt();

//			long st = System.currentTimeMillis();

			preprocessing();
			pw.println(Query());

//			long et = System.currentTimeMillis();
//			System.out.println(et-st);
		}

		pw.close();
	}


	public static void main(String[] args) throws Exception {
		// do not alter this method
		Supermarket ps7 = new Supermarket();
		try{
			ps7.run();
		} catch (Exception e) {
			e.printStackTrace(System.out);
		}
	}
}


class IntegerScanner { // coded by Ian Leow, using any other I/O method is not recommended
	BufferedInputStream bis;
	IntegerScanner(InputStream is) {
		bis = new BufferedInputStream(is, 1000000);
	}

	public int nextInt() {
		int result = 0;
		try {
			int cur = bis.read();
			if (cur == -1)
				return -1;

			while ((cur < 48 || cur > 57) && cur != 45) {
				cur = bis.read();
			}

			boolean negate = false;
			if (cur == 45) {
				negate = true;
				cur = bis.read();
			}

			while (cur >= 48 && cur <= 57) {
				result = result * 10 + (cur - 48);
				cur = bis.read();
			}

			if (negate) {
				return -result;
			}
			return result;
		}
		catch (IOException ioe) {
			return -1;
		}
	}
}


class IntegerPair implements Comparable<IntegerPair> {
	int _first, _second;

	public IntegerPair(int f, int s) {
		_first = f;
		_second = s;
	}

	public int compareTo(IntegerPair o) {
		if (this.first() != o.first())
			return this.first() - o.first();
		else
			return this.second() - o.second();
	}

	int first() { return _first; }
	int second() { return _second; }
}


// inner class Node
class Node {
	int lowerBound;
	int size;

	int[][] matrix;
	byte[] rowIndex; // indicating index of row for matrix in reducedT 
	byte[] colIndex; // indicating index of col for matrix in reducedT 

	int[] rowPenalty;
	int[] colPenalty;
	
	int highestPenaltyRowIndex, highestPenaltyColIndex;
	int maxPenalty;
	
	ArrayList<byte[]> pathHistory = new ArrayList<byte[]>(); // keep track of visited sub-paths

	public Node (int size) {
		this.size = size;

		matrix = new int[size][size];
		rowIndex = new byte[size];
		colIndex = new byte[size];

		rowPenalty = new int[size];
		colPenalty = new int[size];
	}

	/**
	 * 2 * 2 matrix
	 */
	public void calFeasible() {
		if (matrix[0][0] >= 0 && matrix[0][1] >= 0 && matrix[1][0] >= 0 && matrix[1][1] >= 0) {
			lowerBound += Math.min(matrix[0][1] + matrix[1][0], matrix[0][0] + matrix[1][1]);
		} else if (matrix[0][1] >= 0 && matrix[1][0] >= 0) {
			lowerBound = lowerBound + matrix[0][1] + matrix[1][0];
		} else if (matrix[0][0] >= 0 && matrix[1][1] >= 0){
			lowerBound = lowerBound + matrix[0][0] + matrix[1][1];
		}
	}

	/**
	 * Remove entry with highestPenalty to get leftChild
	 */
	public boolean alterForLeftChild() {
		if (this.size == highestPenaltyColIndex) {
			return false;
		}
		
		this.lowerBound += maxPenalty;
		matrix[highestPenaltyRowIndex][highestPenaltyColIndex] = -1;
		
		// reduce row and col value to get update matrix
		for (int i=0; i<size; i++) {
			// update row
			matrix[highestPenaltyRowIndex][i] -= rowPenalty[highestPenaltyRowIndex];
			// update col
			matrix[i][highestPenaltyColIndex] -= colPenalty[highestPenaltyColIndex];
		}
		return true;
	}

	/**
	 * Copy the node input and remove the one row and one col based on input
	 * @param processedNode 
	 * @param highestPenaltyRowIndex2
	 * @param highestPenaltyColIndex2
	 */
	public Node(Node processedNode, int removedRowIndex, int removedColIndex) {
		this.size = processedNode.getSize() - 1;
		this.lowerBound = processedNode.getLowerBound();
		
		matrix = new int[size][size];
		rowIndex = new byte[size];
		colIndex = new byte[size];

		rowPenalty = new int[size];
		colPenalty = new int[size];
		
		int row = 0;
		int col = 0;
		for (int i=0; i<processedNode.getSize(); i++) {
			if (i != removedRowIndex) {
				rowIndex[row] = processedNode.getRowIndex(i);
				col = 0;
				
				for (int j=0; j<processedNode.getSize(); j++) {
					if (j != removedColIndex) {
						matrix[row][col] = processedNode.getMatrix(i, j);
						col++;
					}
				}
				
				row++;
			}
		} // end for loop
		
		col = 0;
		// assign colIndex
		for (int i=0; i<processedNode.getSize(); i++) {
			if (i != removedColIndex) {
				colIndex[col] = processedNode.getColIndex(i);
				col++;
			}
		}
		
		// copy visited sub-paths
		for (byte[] path : processedNode.getPathHistory()) {
			pathHistory.add(path.clone());
		}
		
		// remove unallowed entry by checking for cycle
		removeCycle(processedNode.getRowIndex(removedRowIndex), processedNode.getColIndex(removedColIndex));
	}


	/**
	 * @param row
	 * @param col
	 */
	private void removeCycle(byte row, byte col) {
		// go through row and col, search and mark -ve for entry with (colIndex, rowIndex)
		boolean isFound = false;
		for (int i=0; i<size; i++) {
			if (rowIndex[i] == col) {
				for (int j=0; j<size; j++) {
					if (colIndex[j] == row) {
						matrix[i][j] = -1;
						isFound = true;
					}
				}
			}
		}
		
		if (isFound) {
			byte[] newPath = {row, col};
			pathHistory.add(newPath);
			return;
		} else {
			// trace the pathHistory to find corresponding entry to be removed so that cycle doesn't exists.
			for (int i=0; i<pathHistory.size(); i++) {
				byte[] temp = pathHistory.get(i);
				if (temp[1] == row) {
					temp[1] = col;
					pathHistory.remove(i);
					this.removeCycle(temp[0], temp[1]);
					break;
				} else if (temp[0] == col) {
					temp[0] = row;
					pathHistory.remove(i);
					this.removeCycle(temp[0], temp[1]);
					break;
				}
			}
		}
	}

	/**
	 * Select a node to be branched at based on row and col penalty
	 */
	public void selectNode() {
		this.calRowPenalty();
		this.calColPenalty();
		
		maxPenalty = -1;
		// calculate total penalty for each zero entry and keep track of the index of node with highest penalty
		for (int i=0; i<size; i++) {
			for (int j=0; j<size; j++) {
				if (matrix[i][j] == 0) {
					// check total penalty value
					if (rowPenalty[i] + colPenalty[j] > maxPenalty) {
						// update
						maxPenalty = rowPenalty[i] + colPenalty[j];
						highestPenaltyRowIndex = i;
						highestPenaltyColIndex = j;
					}
				}
			}
		}
	}

	/**
	 * Cal col minimum penalty for each col
	 */
	private void calColPenalty() {
		// col
		for (int i=0; i<size; i++) {
			int numOfZero = 0;
			int minimumPenalty = Integer.MAX_VALUE;
			// row
			for (int j=0; j<size; j++) {
				if (matrix[j][i] == 0) {
					numOfZero++;
					if (numOfZero > 1) {
						minimumPenalty = 0;
						break;
					}
				} else if (matrix[j][i] > 0 && matrix[j][i] < minimumPenalty) {
					minimumPenalty = matrix[j][i];
				}
			}
			colPenalty[i] = minimumPenalty;
		}
	}

	/**
	 * Cal row minimum penalty for each row
	 */
	private void calRowPenalty() {
		// row
		for (int i=0; i<size; i++) {
			int numOfZero = 0;
			int minimumPenalty = Integer.MAX_VALUE;
			// col
			for (int j=0; j<size; j++) {
				if (matrix[i][j] == 0) {
					numOfZero++;
					if (numOfZero > 1) {
						minimumPenalty = 0;
						break;
					}
				} else if (matrix[i][j] > 0 && matrix[i][j] < minimumPenalty) {
					minimumPenalty = matrix[i][j];
				}
			}
			rowPenalty[i] = minimumPenalty;
		}
	}

	/**
	 * subtract col minimal for each col
	 */
	public void reduceCol() {
		for (int i=0; i<size; i++) {
			int colMinimal = finaColMinimal(i);
			// reduce every row entry by rowMinimal
			for (int j=0; j<size; j++) {
				matrix[j][i] -= colMinimal;
			}
			// add rowMinimal to LB
			lowerBound += colMinimal;
		}
	}

	/**
	 * @param i
	 * @return colMinimal
	 */
	private int finaColMinimal(int col) {
		int minimal = Integer.MAX_VALUE;
		for (int i=0; i<size; i++) {
			// find non-negative minimal
			if (matrix[i][col] >= 0 && matrix[i][col] < minimal) {
				minimal = matrix[i][col];
			}
		}
		return minimal;
	}

	/**
	 * subtract row minimal for each row
	 */
	public void reduceRow() {
		for (int i=0; i<size; i++) {
			int rowMinimal = finaRowMinimal(i);
			// reduce every row entry by rowMinimal
			for (int j=0; j<size; j++) {
				matrix[i][j] -= rowMinimal;
			}
			// add rowMinimal to LB
			lowerBound += rowMinimal;
		}
	}

	/**
	 * @param i
	 * @return rowMinimal
	 */
	private int finaRowMinimal(int row) {
		int minimal = Integer.MAX_VALUE;
		for (int i=0; i<size; i++) {
			// find non-negative minimal
			if (matrix[row][i] >= 0 && matrix[row][i] < minimal) {
				minimal = matrix[row][i];
			}
		}
		return minimal;
	}

	/**
	 * @param reducedT
	 */
	public void assignMatrixValue(int[][] originalMatrix) {
		for (int i=0; i<size; i++) {
			for (int j=0; j<size; j++) {
				matrix[i][j] = originalMatrix[i][j];
			}
		}
	}

	public void assignRowIndex(int index, byte value) {
		rowIndex[index] = value;
	}

	public void assignColIndex(int index, byte value) {
		colIndex[index] = value;
	}

	
	public ArrayList<byte[]> getPathHistory() {
		return pathHistory;
	}
	
	public byte getRowIndex(int index) {
		return rowIndex[index];
	}

	public byte getColIndex(int index) {
		return colIndex[index];
	}
	
	public int getSize() {
		return size;
	}
	
	private int getMatrix(int i, int j) {
		return matrix[i][j];
	}

	public int getMaxPenalty() {
		return maxPenalty;
	}

	public int getLowerBound() {
		return lowerBound;
	}

	public int getHighestPenaltyRowIndex() {
		return highestPenaltyRowIndex;
	}

	public int getHighestPenaltyColIndex() {
		return highestPenaltyColIndex;
	}
}