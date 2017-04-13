#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <map>
using namespace std;

void usage(const string prog){
  cout << "usage:" << endl << prog << " <input file>  [outputname]" << endl;
  return;
}

int writeOff(string off, unsigned long long int vC, unsigned long long int fC, vector< vector<float> > sortedVertices, vector< vector<unsigned long long int> > faces)
{
  int res = 1; // function return value
  unsigned long long int i; // loop counter
  
  ofstream offFile(off); // open a new off file to write

  offFile << "OFF" << endl; // first line of an off file
    
  offFile << vC << " " << fC << " " << 0 << endl; // number of vertices + number of faces + number of edges (which could be ignored but at least meshlab requires any number here, so we always write 0)
  for( i = 0 ; i < vC ; i++ ) // loop over the vertices
    {
      vector<float> v = sortedVertices[i]; // extract a vertex
      offFile << v[0] << " " << v[1] << " " << v[2] << endl; // write its components
    }
  for( i = 0 ; i < fC ; i++ ) // loop over the faces
    {
      vector<unsigned long long int> vIndices = faces[i]; // get the 3 indices of the face's vertices
      offFile << "3 " << vIndices[0] << " " << vIndices[1] << " " << vIndices[2] << endl; // write face description: number of vertices (3) + each vertex index
    }

  offFile.close(); // close output file
  if(offFile.is_open())
    cout << "Something went wrong when closing " << off << endl;
  
  res = 0;
  
  return res;
}


int asciiStlToOff(const string stl, const string off){
  
  int res = 1; // function return value
  
  string line; // buffer to read line by line
  string::size_type pos, tmp; // cursor positions to parse a line

  int i, j; // loop counters

  vector<float> n(3); // normal vector
  array< vector<float>, 3 > t; // triangle =  array of 3 vertices
  fill(t.begin(), t.end(), vector<float>(3));

  // we use nested maps for a faster vertex lookup (to merge of duplicated vertices)
  map< float, map< float, map< float, unsigned long long int > > > verticesIndices; // map vertex->index
  // we also use a list of vertices in which their order is known
  vector< vector<float> > sortedVertices; // map index->vertex
  unsigned long long int vC = 0; // vertices counter/index
  bool vIsNew; // true if a vertex is not is the maps above yet

  // a given face is defined by the three indices of its vertices
  vector< vector<unsigned long long int> > faces; // list of faces
  unsigned long long int fC = 0; // faces counter/index

  ifstream stlFile (stl); // stl file to read (input)
  

  if( ! stlFile.is_open() ) cout << "Unable to open file " << stl << endl;
  else{
    getline(stlFile,line); // "solid" line
    getline(stlFile, line); // start reading the first facet (if there is any)
    while( line.compare( 0, 8, "endsolid" ) != 0 ){ // as long as we don't reach the end of the solid description
      fC ++; // we read a new facet (=> first face has index 1 and not 0)
      faces.push_back( vector<unsigned long long int>() );
      pos = line.find("facet") + string("facet normal").length(); // maybe set this after "facet normal" ?
      for( i=0; i<3; i++ ) // parse the normal components
	{
	  n[i] = stof( line.substr(pos), &tmp );
	  pos = pos + tmp;
	}
      getline(stlFile, line); // "outer loop" line
      for( i=0; i<3; i++ ) // parse the 3 vertices
	{
	  getline(stlFile, line); // reads a vertex line
	  pos = line.find("vertex") + string("vertex").length(); // maybe set this after "vertex" ?
	  for( j=0; j<3; j++ ) // parse the vertex components
	    {
	      t[i][j] = stof( line.substr(pos), &tmp );
	      pos = pos + tmp;
	    }
	  float x = t[i][0];
	  float y = t[i][1];
	  float z = t[i][2];
	  if( verticesIndices.count(x) == 0 ) { // x component lookup
	    map< float, map<float, unsigned long long int> > xMap;
	    verticesIndices[x] = xMap;
	    map< float, unsigned long long int>  yMap;
	    verticesIndices[x][y] = yMap;
	    vIsNew = true;
	  }
	  else{
	    if( verticesIndices[x].count(y) == 0 ) { // y component lookup
	      map< float, unsigned long long int>  yMap;
	      verticesIndices[x][y] = yMap;
	      vIsNew = true;
	    }
	    else{	      
	      if( verticesIndices[x][y].count(z) == 0 ) { // z component lookup
		vIsNew = true;
	      }
	      else{
		vIsNew = false;
	      }
	    }
	  }
	  if(vIsNew)
	    {
	      verticesIndices[x][y][z] = vC; // add new vertex index (the current number of vertices, starting at 0)
	      sortedVertices.push_back( vector<float>(t[i]) ); // add new vertex in the list
	      vC++; // increase number of vertices
	    }
	  faces[fC-1].push_back(verticesIndices[x][y][z]); // we add this vertex index to the current face description
	}
      getline(stlFile, line); // "endloop" line
      getline(stlFile, line); // "endfacet" line
      getline(stlFile, line); // start reading either the next facet (if there is any) or the "endsolid" line
    } // end of while loop (last read line is "endsolid"
    
    stlFile.close(); // close the stl file we read
    if(stlFile.is_open())
      cout << "Something went wrong when closing " << stl << endl;
       
    res = writeOff(off, vC, fC, sortedVertices, faces);
  }
  
  return res;
}

int binaryStlToOff(const string stl, const string off){
  
  int res = 1; // function return value
  
  int i; // loop counters

  float *normal = new float [3]; // normal vector
  float *vertex = new float [9]; // triangle = array of 3 vertices = 9 float components
  unsigned short attribute;  // attribute byte count (cf binary stl file format)
  
  float x,y,z; // used as current vector components

  map< float, map< float, map< float, unsigned long long int > > > verticesIndices; // map vertex->index
  vector< vector<float> > sortedVertices; // map index->vertex
  unsigned long long int vC = 0; // vertices counter/index
  bool vIsNew; // true if a vertex is not is the maps above

  vector< vector<unsigned long long int> > faces; // map index->faxe
  unsigned long long int fC = 0; // faces counter/index
  unsigned int nTriangles; // number of triangles written in binary STL files, should be equals to fC in the end;
  
  ifstream sf (stl, ifstream::binary); // stl file to read
  streampos begin, end;

  if( ! sf.is_open() ) cout << "Unable to open file " << stl << endl;
  else{

    begin = sf.tellg();
    sf.seekg(0,ios::end);
    end = sf.tellg();
    unsigned int totalSize = end - begin;
    unsigned int headerSize = 80 + 4; // 80-byte header and 4-byte number of triangles
    unsigned int allBlocksSize = totalSize - headerSize; // blocks representing faces = the rest of the file
    unsigned int blockSize = 4*12+2; // size of one block: 4*3 4-byte floats (normal vector and 3 vertices) + 2-byte attribute count
    if( allBlocksSize % blockSize !=0 ) {
      cout << "Inconsistent file size: face blocks count should be a whole number ! Aborting..." << endl;
      return 1;
    }
    unsigned int nFaces = allBlocksSize / blockSize;
    sf.seekg(0,ios::beg);
    char *header = new char [80]; // To read and store the header
    sf.read( header, 80 ); // We skip the 80-byte header
    sf.read( (char*) &nTriangles, 4 ); // We read the UINT32 (4 bytes) number of triangles
    if( nTriangles != nFaces ) {
      cout << "Warning: Inconsistent number of triangles read ! The file says " << nTriangles << " triangles, but has " << nFaces << " face blocks."  << endl;
      return 1;
    }
    for( fC = 0; fC < nFaces; fC++ ) { // for each new facet block
      sf.read( (char*) normal, 12 ); // We read the normal vector: 3 REAL32 (float, 4 bytes) components
      sf.read( (char*) vertex, 12 ); // We read the x vertex: 3 REAL32 (float, 4 bytes) components 
      sf.read( (char*) &(vertex[3]), 12 ); // We read the y vertex: 3 REAL32 (float, 4 bytes) components 
      sf.read( (char*) &(vertex[6]), 12 ); // We read the z vertex: 3 REAL32 (float, 4 bytes) components 
      sf.read( (char*) &attribute, 2 ); // We read the UINT16 (unsigned short, 2 bytes) attribute byte count
      faces.push_back( vector<unsigned long long int>() ); // We create a new (empty) face
      for( i=0; i<3; i++ ) {// parse each vertex components
	x =  vertex[3*i];
	y =  vertex[3*i+1];
	z =  vertex[3*i+2];
	if( verticesIndices.count(x) == 0 ) { // x component lookup
	  map< float, map<float, unsigned long long int> > xMap; // create a new map for this x value
	  verticesIndices[x] = xMap;
	  map< float, unsigned long long int>  yMap; // create a new map for this y value
	  verticesIndices[x][y] = yMap;
	  vIsNew = true;
	}
	else{
	  if( verticesIndices[x].count(y) == 0 ) { // y component lookup
	    map< float, unsigned long long int>  yMap; // create a new map for this y value
	    verticesIndices[x][y] = yMap;
	    vIsNew = true;
	  }
	  else{	      
	    if( verticesIndices[x][y].count(z) == 0 ) { // z component lookup
	      vIsNew = true;
	    }
	    else{
	      vIsNew = false;
	    }
	  }
	}
	if(vIsNew) {
	  verticesIndices[x][y][z] = vC; // add new vertex index (map z to vC)
	  sortedVertices.push_back( vector<float>( &(vertex[3*i]), &(vertex[3*i]) + 3 ) ); // add new vertex to the list
	  vC++; // increase vertex count
	}
	faces[fC].push_back(verticesIndices[x][y][z]); // add the vertex index to the current face description	
      }
    }

    sf.close(); // close the stl file we read
    if(sf.is_open())
      cout << "Something went wrong when closing " << stl << endl;
    res = writeOff(off, vC, fC, sortedVertices, faces);
  }

  return res;
}


int main (int argc, char **argv) {

  if( argc != 2 && argc != 3 ){ // check argument count (should be 2 or 3)
    string prog(argv[0]);
    usage(string(prog));
    return 0;
  }

  string stl(argv[1]);
  string off; // output name
  if(argc == 2) off = stl.substr( 0, stl.rfind('.') ) + ".off"; // if no output name was supplied: use input name + 'off' extension
  if(argc == 3) off = string(argv[2]); 
  string line; // buffer to read line by line

  int res = 1; // main return value
  
  ifstream  stlFile (stl);

  if( ! stlFile.is_open() ) cout << "Unable to open file " << stl << endl;
  else{
    getline(stlFile,line);
    stlFile.close();
    if( line.compare( 0, 5, "solid" ) == 0 ){ // if the file is ASCII defined
      cout << "Input: "<< stl << " (ASCII STL) - Output: " << off << endl;
      res = asciiStlToOff(stl, off); // call ASCII STL to OFF converter
    }
    else{ // else
      cout << "Input: " << stl << ": (binary STL) - Output: " << off << endl;
      res = binaryStlToOff(stl, off); // process it as a binary STL.
    }
  }
  return res;
}
