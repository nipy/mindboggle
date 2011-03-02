#include "SurfaceIO.h"
#include "pathfinder.h"

#define LENGTH_THRESHOLD 15  // A LOWER threshold of Sulci Length
#define MAX_SULCI_LENGTH 16384 // the maximum length allowed 
#define EPSILON 0.0001


void RecordSegment(int v0, int v1, int v2,
			float emax0, float emax1, float emax2,
			float kmax0, float kmax1, float kmax2,
			float thresh, bool to_center, Segment *segment, Surface *surface, int faceIndex, bool isStrictSegment)
{
	int edgeS, edgeE;
	float w10 = fabs(emax0) / (fabs(emax0) + fabs(emax1));
	float w01 = 1.0f - w10;
	float k01 = fabs(w01 * kmax0 + w10 * kmax1);

	Fvector3d p01, p12;
	Fvector3d pdir01, pdir12;
	float k12;

	p01 =  Fvector3dADDFvector3d( ScalarMultiplyFvector3d(w01,surface->vertex[v0]),ScalarMultiplyFvector3d(w10,surface->vertex[v1]));
	
	pdir01 =  Fvector3dADDFvector3d( ScalarMultiplyFvector3d(w01,surface->pdir2[v0]),ScalarMultiplyFvector3d(w10,surface->pdir2[v1]));
	
	FindEdgeID(surface, v0, v1, &edgeS);

	if (to_center) 
	{
		// Connect first point to center of triangle
		p12.x = (surface->vertex[v0].x + surface->vertex[v1].x + surface->vertex[v2].x) / 3.0f;
		p12.y = (surface->vertex[v0].y + surface->vertex[v1].y + surface->vertex[v2].y) / 3.0f;
		p12.z = (surface->vertex[v0].z + surface->vertex[v1].z + surface->vertex[v2].z) / 3.0f;

		pdir12.x = (surface->pdir2[v0].x + surface->pdir2[v1].x + surface->pdir2[v2].x) / 3.0f;
		pdir12.y = (surface->pdir2[v0].y + surface->pdir2[v1].y + surface->pdir2[v2].y) / 3.0f;
		pdir12.z = (surface->pdir2[v0].z + surface->pdir2[v1].z + surface->pdir2[v2].z) / 3.0f;

		k12 = fabs(kmax0 + kmax1 + kmax2) / 3.0f;

		edgeE = -2;
	} 
	else 
	{
		// Connect first point to second one (on next edge)
		float w21 = fabs(emax1) / (fabs(emax1) + fabs(emax2));
		float w12 = 1.0f - w21;
		
		p12 = Fvector3dADDFvector3d( ScalarMultiplyFvector3d(w12,surface->vertex[v1]),ScalarMultiplyFvector3d(w21,surface->vertex[v2]));
		pdir12 = Fvector3dADDFvector3d( ScalarMultiplyFvector3d(w12,surface->pdir2[v1]),ScalarMultiplyFvector3d(w21,surface->pdir2[v2]));

		k12 = fabs(w12 * kmax1 + w21 * kmax2);

		FindEdgeID(surface, v1, v2, &edgeE);
	}

	if(k01 < thresh || k12 < thresh)
		return;

	NormalizeFvector3d(&pdir01);
	NormalizeFvector3d(&pdir12);
	
	segment->isStrictSegment.push_back(isStrictSegment); 
	segment->start.push_back(p01);
	segment->end.push_back(p12);
	segment->pdirS.push_back(pdir01);
	segment->pdirE.push_back(pdir12);
	segment->curvS.push_back(k01);
	segment->curvE.push_back(k12);
	segment->edgeS.push_back(edgeS);
	segment->edgeE.push_back(edgeE);
	segment->face.push_back(faceIndex);
	segment->Face2Segment[faceIndex].push_back(segment->n);
	segment->endIsCenter.push_back(to_center);

	segment->n++;

	return;
}


void ComputeSulciSegment(Surface *surface, Params *params, Segment *segment, bool isGyri)
{
	int i;
	float thresh = 0.1;
	float rv_sign = isGyri ? 1.0f : -1.0f;

	segment->n = 0;

	segment->Face2Segment.resize(surface->faceNum);

	for(i=0; i<surface->faceNum; i++)
	{
		int v0, v1, v2;

		v0 = surface->faces[i].v[0];
		v1 = surface->faces[i].v[1];
		v2 = surface->faces[i].v[2];

		if(isGyri)
		{
			if ((surface->curv1[v0] <= 0.0f) || (surface->curv1[v1] <= 0.0f) || (surface->curv1[v2] <= 0.0f))
				continue;
		}
		else
		{
			if ((surface->curv1[v0] >= 0.0f) || (surface->curv1[v1] >= 0.0f) || (surface->curv1[v2] >= 0.0f))
				continue;
		}

		float emax0, emax1, emax2;
		
		emax0 =  surface->dcurv[0][v0];
		emax1 =  surface->dcurv[0][v1];
		emax2 =  surface->dcurv[0][v2];

		Fvector3d tmax0, tmax1, tmax2;
		
		tmax0 =  ScalarMultiplyFvector3d(rv_sign*surface->dcurv[0][v0], surface->pdir1[v0]);
		tmax1 =  ScalarMultiplyFvector3d(rv_sign*surface->dcurv[0][v1], surface->pdir1[v1]);
		tmax2 =  ScalarMultiplyFvector3d(rv_sign*surface->dcurv[0][v2], surface->pdir1[v2]);

		bool z01, z12, z20; // used to select candidate sulci segment
		bool zz01, zz12, zz20; // used to select staring search sulci segment
		
		z01 = (Fvector3dDOTFvector3d(tmax0, tmax1) <= 0.0f);
		z12 = (Fvector3dDOTFvector3d(tmax1, tmax2) <= 0.0f);
		z20 = (Fvector3dDOTFvector3d(tmax2, tmax0) <= 0.0f);

		if (z01 + z12 + z20 >= 2)
		{
			Fvector3d p0, p1, p2;
			
			p0.x = surface->vertex[v0].x;
			p0.y = surface->vertex[v0].y;
			p0.z = surface->vertex[v0].z;

			p1.x = surface->vertex[v1].x;
			p1.y = surface->vertex[v1].y;
			p1.z = surface->vertex[v1].z;

			p2.x = surface->vertex[v2].x;
			p2.y = surface->vertex[v2].y;
			p2.z = surface->vertex[v2].z;

			// Check whether we have the correct flavor of extremum:
			// Is the curvature increasing along the edge?
			zz01 = z01 &&  (Fvector3dDOTFvector3d(tmax0, Fvector3dMINUSFvector3d(p1,p0)) >= 0.0f || 
				Fvector3dDOTFvector3d(tmax1, Fvector3dMINUSFvector3d(p1,p0)) <= 0.0f);

			zz12 = z12 &&  (Fvector3dDOTFvector3d(tmax1, Fvector3dMINUSFvector3d(p2,p1)) >= 0.0f || 
				Fvector3dDOTFvector3d(tmax2, Fvector3dMINUSFvector3d(p2,p1)) <= 0.0f);

			zz20 = z20 &&  (Fvector3dDOTFvector3d(tmax2, Fvector3dMINUSFvector3d(p0,p2)) >= 0.0f || 
				Fvector3dDOTFvector3d(tmax0, Fvector3dMINUSFvector3d(p0,p2)) <= 0.0f);

			// strat segment
			if (zz01 + zz12 + zz20 >= 2)
			{
				float kmax0 = surface->curv1[v0];
				float kmax1 = surface->curv1[v1];
				float kmax2 = surface->curv1[v2];

				int faceIndex = i;

				if (!zz01) 
				{
					RecordSegment(v1, v2, v0, 
						emax1, emax2, emax0, 
						kmax1, kmax2, kmax0, 
						thresh, false, segment, surface, faceIndex, true);
				}
				else if (!zz12) 
				{
					RecordSegment(v2, v0, v1,
						emax2, emax0, emax1,
						kmax2, kmax0, kmax1,
						thresh, false, segment, surface, faceIndex, true);
				} 
				else if (!zz20) 
				{
					RecordSegment(v0, v1, v2,
						emax0, emax1, emax2,
						kmax0, kmax1, kmax2,
						thresh, false, segment, surface, faceIndex, true);
				} 
				else 
				{
					// All three edges have crossings -- connect all to center
					RecordSegment(v1, v2, v0,
						emax1, emax2, emax0,
						kmax1, kmax2, kmax0,
						thresh, true, segment, surface, faceIndex, true);
					
					RecordSegment(v2, v0, v1,
						emax2, emax0, emax1,
						kmax2, kmax0, kmax1,
						thresh, true, segment, surface, faceIndex, true);
					
					RecordSegment(v0, v1, v2,
						emax0, emax1, emax2,
						kmax0, kmax1, kmax2,
						thresh, true, segment, surface, faceIndex, true);
				}
			}
			//candidate segment
			else
			{
				float kmax0 = surface->curv1[v0];
				float kmax1 = surface->curv1[v1];
				float kmax2 = surface->curv1[v2];

				int faceIndex = i;

				if (!z01) 
				{
					RecordSegment(v1, v2, v0, 
						emax1, emax2, emax0, 
						kmax1, kmax2, kmax0, 
						thresh, false, segment, surface, faceIndex, false);
				}
				else if (!z12) 
				{
					RecordSegment(v2, v0, v1,
						emax2, emax0, emax1,
						kmax2, kmax0, kmax1,
						thresh, false, segment, surface, faceIndex, false);
				} 
				else if (!z20) 
				{
					RecordSegment(v0, v1, v2,
						emax0, emax1, emax2,
						kmax0, kmax1, kmax2,
						thresh, false, segment, surface, faceIndex, false);
				} 
				else 
				{
					// All three edges have crossings -- connect all to center
					RecordSegment(v1, v2, v0,
						emax1, emax2, emax0,
						kmax1, kmax2, kmax0,
						thresh, true, segment, surface, faceIndex, false);
					
					RecordSegment(v2, v0, v1,
						emax2, emax0, emax1,
						kmax2, kmax0, kmax1,
						thresh, true, segment, surface, faceIndex, false);
					
					RecordSegment(v0, v1, v2,
						emax0, emax1, emax2,
						kmax0, kmax1, kmax2,
						thresh, true, segment, surface, faceIndex, false);
				}
			}
		}
	}

	return;
}


void WriteSegmentVTK(char *fname, Segment *segment)
{
	FILE *fp;

	int i, t, pointNum;

	fp = fopen(fname,"wt");

	fprintf(fp,"# vtk DataFile Version 3.0\n");
	fprintf(fp,"Segment\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp,"POINTS %d float\n",2*segment->n);
	
	for(i=0; i<segment->n; i++)
	{
		fprintf(fp,"%f %f %f\n",segment->start[i].x, segment->start[i].y, segment->start[i].z);
		fprintf(fp,"%f %f %f\n",segment->end[i].x, segment->end[i].y, segment->end[i].z);
	}

	fprintf(fp,"CELLS %d %d\n", segment->n, segment->n*3);
	
	for(i=0; i<segment->n; i++)
	{
		fprintf(fp,"2 %d %d\n",i*2, i*2+1);
	}

	fprintf(fp,"CELL_TYPES %d\n", segment->n);
	for(i=0; i<segment->n; i++)
	{
		fprintf(fp,"4\n");
	}

	fprintf(fp,"POINT_DATA %d\n",segment->n*2);
	fprintf(fp,"SCALARS curv1 float\n");
	fprintf(fp,"LOOKUP_TABLE curv1Table\n");
	
	for(i=0; i<segment->n; i++)
	{
		fprintf(fp,"%f\n",segment->curvS[i]);
		fprintf(fp,"%f\n",segment->curvE[i]);
	}

	fprintf(fp,"SCALARS isStartSegment int\n");
	fprintf(fp,"LOOKUP_TABLE isStartSegmentTable\n");
	
	for(i=0; i<segment->n; i++)
	{
		if(segment->isStartSegment[i] > 0)
		{
			fprintf(fp,"255\n");
			fprintf(fp,"255\n");
		}
		else
		{
			fprintf(fp,"0\n");
			fprintf(fp,"0\n");
		}
	}

	fprintf(fp,"SCALARS isStrictSegment int\n");
	fprintf(fp,"LOOKUP_TABLE isStrictSegmentTable\n");
	
	for(i=0; i<segment->n; i++)
	{
		if(segment->isStrictSegment[i] > 0)
		{
			fprintf(fp,"255\n");
			fprintf(fp,"255\n");
		}
		else
		{
			fprintf(fp,"0\n");
			fprintf(fp,"0\n");
		}
	}
	
	fprintf(fp,"VECTORS pdir1 float\n");
	for(i=0; i<segment->n; i++)
	{
		fprintf(fp,"%f %f %f\n",segment->pdirS[i].x, segment->pdirS[i].y, segment->pdirS[i].z);
		fprintf(fp,"%f %f %f\n",segment->pdirE[i].x, segment->pdirE[i].y, segment->pdirE[i].z);
	}

	fclose(fp);

	return;
}


bool IsEqual(Fvector3d vector1, Fvector3d vector2)
{
	float diff;

	diff = sqrt( (vector1.x - vector2.x)*(vector1.x - vector2.x) 
		+ (vector1.y - vector2.y)*(vector1.y - vector2.y)
		+ (vector1.z - vector2.z)*(vector1.z - vector2.z));

	return (diff < EPSILON);
}

void FindAdjacentSegmentPairs(Surface *surface, Segment *segment)
{
	int i, j, processed;
	int currentSegment, nextSegment;
	int currentFace, nextFace;
	
	segment->adjacentSegment.resize(segment->n);

	for(processed=0; processed<segment->n; processed++)
	{
		currentSegment = processed;
		currentFace = segment->face[processed];

		// 3 acorss_edge faces
		for(i=0; i<3; i++)
		{
			nextFace = surface->across_edge[currentFace].v[i];

			if( nextFace >0 )
			{
				for(j=0; j<segment->Face2Segment[nextFace].size(); j++)
				{
					nextSegment = segment->Face2Segment[nextFace][j];

					segment->adjacentSegment[currentSegment].push_back(nextSegment);
				}
			}
		}

		int count = 0;

		for(i=0; i<segment->adjacentSegment[currentSegment].size(); i++)
		{
			if( IsEqual(segment->start[currentSegment],segment->start[segment->adjacentSegment[currentSegment][i]])
				|| IsEqual(segment->start[currentSegment],segment->end[segment->adjacentSegment[currentSegment][i]])
				|| IsEqual(segment->end[currentSegment],segment->start[segment->adjacentSegment[currentSegment][i]])
				|| IsEqual(segment->end[currentSegment],segment->end[segment->adjacentSegment[currentSegment][i]]) )
				count++;
		}

		// startimg segment 
		if(count == 1 && segment->Face2Segment[currentFace].size() == 1)
			segment->isStartSegment.push_back(1);
		// single segment
		else if(count == 0)
			segment->isStartSegment.push_back(2);
		// connecting segment
		else
			segment->isStartSegment.push_back(0);
	}
	
	return;
}


void SwitchSegmentStartEnd(Segment *segment, int i)
{
	Fvector3d vector;
	float value;
	int valueInt; 

	value = segment->curvS[i];
	segment->curvS[i] = segment->curvE[i];
	segment->curvE[i] = value;

	valueInt = segment->edgeS[i];
	segment->edgeS[i] = segment->edgeE[i];
	segment->edgeE[i] = valueInt;

	vector.x = segment->start[i].x;
	vector.y = segment->start[i].y;
	vector.z = segment->start[i].z;

	segment->start[i].x = segment->end[i].x;
	segment->start[i].y = segment->end[i].y;
	segment->start[i].z = segment->end[i].z;

	segment->end[i].x = vector.x;
	segment->end[i].y = vector.y;
	segment->end[i].z = vector.z;
    
	vector.x = segment->pdirS[i].x;
	vector.y = segment->pdirS[i].y;
	vector.z = segment->pdirS[i].z;

	segment->pdirS[i].x = segment->pdirE[i].x;
	segment->pdirS[i].y = segment->pdirE[i].y;
	segment->pdirS[i].z = segment->pdirE[i].z;

	segment->pdirE[i].x = vector.x;
	segment->pdirE[i].y = vector.y;
	segment->pdirE[i].z = vector.z;

	return;
}

void ConnectSulciFromSegment(SulciTrackingOut *out, Segment *segment)
{
	Surface *surface; // a pointer defined for convenience, points to the surface storage space 
	Sulci contour;
	unsigned char *segmentFlag; // used to label whether a segment is processed
	
	int processed;	
	int direction;
	int curveLength;
	int currentSegment, nextSegment;
	bool selectStart, selectEnd;
	bool existingStrictSegment;
	int i;

	// a pointer to the surface 
	surface = &(out->surface);
	
	out->numberOfSulci = 0;
	out->sulci = (Sulci*)malloc(sizeof(Sulci)*(out->numberOfSulci+1));
	
	segmentFlag = (unsigned char *)malloc(sizeof(unsigned char)*segment->n);
	memset(segmentFlag,0,sizeof(unsigned char)*segment->n);

	contour.vertex = (Fvector3d *)malloc(sizeof(Fvector3d) * MAX_SULCI_LENGTH);
	contour.curv1 = (float *)malloc(sizeof(float) * MAX_SULCI_LENGTH);
	contour.pdir2 = (Fvector3d *)malloc(sizeof(Fvector3d) * MAX_SULCI_LENGTH);
	contour.segmentID =  (int *)malloc(sizeof(int) * MAX_SULCI_LENGTH);
	contour.faceID =  (int *)malloc(sizeof(int) * MAX_SULCI_LENGTH);
	contour.edgeID = (int *)malloc(sizeof(int) * MAX_SULCI_LENGTH);
	contour.isCenter = (bool *)malloc(sizeof(bool) * MAX_SULCI_LENGTH);

	FindAdjacentSegmentPairs(surface, segment);

	// loop while there remain unprocessed segment
	for(processed=0; processed<segment->n; processed++)
	{
		contour.n = 0;
		contour.junction = 0;
		existingStrictSegment = false;
		
		if( segment->isStartSegment[processed] == 1 )
		{
			// skip this point if it has already been processed 
			if(segmentFlag[processed])
				continue;
		
			// set the flag to 0 so that it is not processed again 
			segmentFlag[processed] = 1;
		
			// Initialise the contour 
			curveLength = 0;
	
			if( IsEqual(segment->start[processed], segment->start[segment->adjacentSegment[processed][0]]) 
			||IsEqual(segment->start[processed], segment->end[segment->adjacentSegment[processed][0]]) )
			{
				SwitchSegmentStartEnd(segment, processed);

				if( IsEqual(segment->end[processed], segment->end[segment->adjacentSegment[processed][0]]) )
				{
					SwitchSegmentStartEnd(segment, segment->adjacentSegment[processed][0]);
				}

				contour.curv1[curveLength] = segment->curvS[processed];

				contour.pdir2[curveLength].x = segment->pdirS[processed].x;
				contour.pdir2[curveLength].y = segment->pdirS[processed].y;
				contour.pdir2[curveLength].z = segment->pdirS[processed].z;

				contour.vertex[curveLength].x = segment->start[processed].x;
				contour.vertex[curveLength].y = segment->start[processed].y;
				contour.vertex[curveLength].z = segment->start[processed].z;
		
				contour.segmentID[curveLength] = processed;
				contour.faceID[curveLength] = segment->face[processed];
				contour.edgeID[curveLength] = segment->edgeS[processed];
				contour.isCenter[curveLength] = segment->endIsCenter[processed];

				curveLength++;

				contour.curv1[curveLength] = segment->curvE[processed];

				contour.pdir2[curveLength].x = segment->pdirE[processed].x;
				contour.pdir2[curveLength].y = segment->pdirE[processed].y;
				contour.pdir2[curveLength].z = segment->pdirE[processed].z;

				contour.vertex[curveLength].x = segment->end[processed].x;
				contour.vertex[curveLength].y = segment->end[processed].y;
				contour.vertex[curveLength].z = segment->end[processed].z;
					
				contour.segmentID[curveLength] = processed;
				contour.faceID[curveLength] = segment->face[processed];
				contour.edgeID[curveLength] = segment->edgeE[processed];
				contour.isCenter[curveLength] = segment->endIsCenter[processed];

				curveLength++;

				if(segment->isStrictSegment[processed])
					existingStrictSegment = true;

				currentSegment = segment->adjacentSegment[processed][0];
			}
			else if( IsEqual(segment->end[processed], segment->start[segment->adjacentSegment[processed][0]]) 
			||IsEqual(segment->end[processed], segment->end[segment->adjacentSegment[processed][0]]) )
			{
				
				if( IsEqual(segment->end[processed], segment->end[segment->adjacentSegment[processed][0]]) )
				{
					SwitchSegmentStartEnd(segment, segment->adjacentSegment[processed][0]);
				}
			
				contour.curv1[curveLength] = segment->curvS[processed];

				contour.pdir2[curveLength].x = segment->pdirS[processed].x;
				contour.pdir2[curveLength].y = segment->pdirS[processed].y;
				contour.pdir2[curveLength].z = segment->pdirS[processed].z;

				contour.vertex[curveLength].x = segment->start[processed].x;
				contour.vertex[curveLength].y = segment->start[processed].y;
				contour.vertex[curveLength].z = segment->start[processed].z;
				
				contour.segmentID[curveLength] = processed;
				contour.faceID[curveLength] = segment->face[processed];
				contour.edgeID[curveLength] = segment->edgeS[processed];
				contour.isCenter[curveLength] = segment->endIsCenter[processed];

				curveLength++;

				contour.curv1[curveLength] = segment->curvE[processed];

				contour.pdir2[curveLength].x = segment->pdirE[processed].x;
				contour.pdir2[curveLength].y = segment->pdirE[processed].y;
				contour.pdir2[curveLength].z = segment->pdirE[processed].z;

				contour.vertex[curveLength].x = segment->end[processed].x;
				contour.vertex[curveLength].y = segment->end[processed].y;
				contour.vertex[curveLength].z = segment->end[processed].z;
				
				contour.segmentID[curveLength] = processed;
				contour.faceID[curveLength] = segment->face[processed];
				contour.edgeID[curveLength] = segment->edgeE[processed];
				contour.isCenter[curveLength] = segment->endIsCenter[processed];
				
				curveLength++;

				if(segment->isStrictSegment[processed])
					existingStrictSegment = true;

				currentSegment = segment->adjacentSegment[processed][0];
			}
			else
			{
				contour.n = 0;
				continue;
			}

			int flag=1;
			
			contour.junction = 0;

			while(flag)
			{      
				segmentFlag[currentSegment] = 1;
				{
					contour.curv1[curveLength] = segment->curvE[currentSegment ];
				
					contour.vertex[curveLength].x = segment->end[currentSegment ].x;
					contour.vertex[curveLength].y = segment->end[currentSegment ].y;
					contour.vertex[curveLength].z = segment->end[currentSegment ].z;
				
					contour.pdir2[curveLength].x = segment->pdirE[currentSegment ].x;
					contour.pdir2[curveLength].y = segment->pdirE[currentSegment ].y;
					contour.pdir2[curveLength].z = segment->pdirE[currentSegment ].z;

					if(!existingStrictSegment && segment->isStrictSegment[currentSegment ])
						existingStrictSegment = true;
					
					contour.segmentID[curveLength] = currentSegment ;
					contour.faceID[curveLength] = segment->face[currentSegment ];
					contour.edgeID[curveLength] = segment->edgeE[currentSegment ];
					contour.isCenter[curveLength] = segment->endIsCenter[currentSegment ];

					curveLength++;
				}
				
				flag = 0;

				for(int iter=0; iter<segment->adjacentSegment[currentSegment].size(); iter++)
				{
					if(!segmentFlag[segment->adjacentSegment[currentSegment][iter]])
					{
						if( IsEqual(segment->end[currentSegment],segment->end[segment->adjacentSegment[currentSegment][iter]]) )
						{
							SwitchSegmentStartEnd(segment, segment->adjacentSegment[currentSegment][iter]);
							nextSegment = segment->adjacentSegment[currentSegment][iter];
							flag = 1;
							break;
						}
						else if( IsEqual(segment->end[currentSegment],segment->start[segment->adjacentSegment[currentSegment][iter]]) )
						{
							nextSegment = segment->adjacentSegment[currentSegment][iter];
							flag = 1;
							break;
						}
					}
				}
				if(flag)
					currentSegment = nextSegment;
			}
			
			// end of the Contour reached 
			contour.n = curveLength;
		}
		
		// check if the length of the sulcus is greater than the threshhold
		if( (existingStrictSegment && contour.n >= LENGTH_THRESHOLD) || contour.junction == 1) 
		{
			// record the curve
			if(out->sulci)
				out->sulci = (Sulci*)realloc(out->sulci,sizeof(Sulci)*(out->numberOfSulci+1));
			else
				out->sulci = (Sulci*)malloc(sizeof(Sulci)*(out->numberOfSulci+1));
			
			// allocate memory
			out->sulci[out->numberOfSulci].n = contour.n;
			out->sulci[out->numberOfSulci].curv1 = (float *)malloc(sizeof(float)*out->sulci[out->numberOfSulci].n);
			out->sulci[out->numberOfSulci].vertex = (Fvector3d *)malloc(sizeof(Fvector3d)*out->sulci[out->numberOfSulci].n);
			out->sulci[out->numberOfSulci].pdir2 = (Fvector3d *)malloc(sizeof(Fvector3d)*out->sulci[out->numberOfSulci].n);
			out->sulci[out->numberOfSulci].segmentID = (int *)malloc(sizeof(int)*out->sulci[out->numberOfSulci].n);
			out->sulci[out->numberOfSulci].faceID = (int *)malloc(sizeof(int)*out->sulci[out->numberOfSulci].n);
			out->sulci[out->numberOfSulci].vertexID = (int *)malloc(sizeof(int)*out->sulci[out->numberOfSulci].n);
			out->sulci[out->numberOfSulci].edgeID = (int *)malloc(sizeof(int)*out->sulci[out->numberOfSulci].n);
			out->sulci[out->numberOfSulci].isJunction = (bool *)malloc(sizeof(bool)*out->sulci[out->numberOfSulci].n);
			out->sulci[out->numberOfSulci].isCenter = (bool *)malloc(sizeof(bool)*out->sulci[out->numberOfSulci].n);

			float len = 0.0;
			for(i=0;i<contour.n;i++)
			{
				out->sulci[out->numberOfSulci].curv1[i] = contour.curv1[i];

				out->sulci[out->numberOfSulci].vertex[i].x = contour.vertex[i].x;
				out->sulci[out->numberOfSulci].vertex[i].y = contour.vertex[i].y;
				out->sulci[out->numberOfSulci].vertex[i].z = contour.vertex[i].z;

				out->sulci[out->numberOfSulci].pdir2[i].x = contour.pdir2[i].x;
				out->sulci[out->numberOfSulci].pdir2[i].y = contour.pdir2[i].y;
				out->sulci[out->numberOfSulci].pdir2[i].z = contour.pdir2[i].z;

				out->sulci[out->numberOfSulci].segmentID[i] = contour.segmentID[i];
				out->sulci[out->numberOfSulci].faceID[i] = contour.faceID[i];
				out->sulci[out->numberOfSulci].edgeID[i] = contour.edgeID[i];
				out->sulci[out->numberOfSulci].vertexID[i] = -1;
				out->sulci[out->numberOfSulci].isCenter[i] = contour.isCenter[i];
				out->sulci[out->numberOfSulci].isJunction[i] = false;

				if(i > 0)
					len += DistanceBetweenTwoVertices(out->sulci[out->numberOfSulci].vertex[i-1], out->sulci[out->numberOfSulci].vertex[i]);
			}
			out->sulci[out->numberOfSulci].color = out->numberOfSulci;
			out->sulci[out->numberOfSulci].junction = 0;
			out->sulci[out->numberOfSulci].len = len;

			// increment the number of curves 
			out->numberOfSulci++;
		}		

	}

	printf("Total %d sulci\n",out->numberOfSulci);

	for(i=0; i<out->numberOfSulci; i++)
	{
		printf("Sulci %d: length %d\n",i,out->sulci[i].n);
	}

	free(segmentFlag);
	free(contour.curv1);
	free(contour.vertex);
	free(contour.pdir2);
	free(contour.faceID);
	free(contour.edgeID);
	free(contour.segmentID);
	free(contour.isCenter);

	return;
}


void RemoveCandiateStartingSegment(SulciTrackingOut *out, Segment *segment)
{
	int i, j, k;
	int count = 0;

	for(i=0; i<out->numberOfSulci; i++)
	{
		printf("Sulci %d Length: %d\n",i, out->sulci[i].n);

		j = 0;

		while(!segment->isStrictSegment[out->sulci[i].segmentID[j]])
			j++;

		if(j > 0)
		{
			for(k=j; k<out->sulci[i].n; k++)
			{
				out->sulci[i].curv1[k-j] = out->sulci[i].curv1[k];

				out->sulci[i].vertex[k-j].x = out->sulci[i].vertex[k].x;
				out->sulci[i].vertex[k-j].y = out->sulci[i].vertex[k].y;
				out->sulci[i].vertex[k-j].z = out->sulci[i].vertex[k].z;

				out->sulci[i].pdir2[k-j].x = out->sulci[i].pdir2[k].x;
				out->sulci[i].pdir2[k-j].y = out->sulci[i].pdir2[k].y;
				out->sulci[i].pdir2[k-j].z = out->sulci[i].pdir2[k].z;

				out->sulci[i].segmentID[k-j] = out->sulci[i].segmentID[k];
				out->sulci[i].faceID[k-j] = out->sulci[i].faceID[k];
				out->sulci[i].edgeID[k-j] = out->sulci[i].edgeID[k];
				out->sulci[i].vertexID[k-j] = out->sulci[i].vertexID[k];
				out->sulci[i].isCenter[k-j] = out->sulci[i].isCenter[k];
				out->sulci[i].isJunction[k-j] = out->sulci[i].isJunction[k];
			}
		}
		
		out->sulci[i].n = out->sulci[i].n - j;

		j = out->sulci[i].n-1;

		while(!segment->isStrictSegment[out->sulci[i].segmentID[j]])
			j--;

		out->sulci[i].n = j+1;

		printf("Sulci %d Length: %d\n",i, out->sulci[i].n);
	}

	return;
}

void InitializeColor2SulciTable(SulciTrackingOut *out)
{
	int i, j, k;
	int count;
	
	out->color2Sulci.resize(out->numberOfSulci);

	for(i=0; i<out->numberOfSulci; i++)
		out->color2Sulci[out->sulci[i].color].push_back(i);
	
	for(i=0; i<out->numberOfSulci; i++)
	for(j=i+1; j<out->numberOfSulci; j++)
	{
		if(out->sulci[i].faceID[out->sulci[i].n-1] == out->sulci[j].faceID[out->sulci[j].n-1] && out->sulci[i].color != out->sulci[j].color)
		{
			if(out->color2Sulci[out->sulci[i].color].size() >= out->color2Sulci[out->sulci[j].color].size())
			{
				count = out->color2Sulci[out->sulci[j].color].size();
				
				for(k=0; k<count; k++)
				{
					out->color2Sulci[out->sulci[i].color].push_back(out->color2Sulci[out->sulci[j].color][k]);
					out->color2Sulci[out->sulci[j].color].pop_back();
					out->sulci[out->color2Sulci[out->sulci[j].color][k]].color = out->sulci[i].color;
				}
			}
			else
			{
				count = out->color2Sulci[out->sulci[i].color].size();

				for(k=0; k<count; k++)
				{
					out->color2Sulci[out->sulci[j].color].push_back(out->color2Sulci[out->sulci[i].color][k]);
					out->color2Sulci[out->sulci[i].color].pop_back();
					out->sulci[out->color2Sulci[out->sulci[i].color][k]].color = out->sulci[j].color;					
				}
			}
		}
	}

	return;
}


void RemoveSingleCenterPoint(SulciTrackingOut *out)
{
	int i, j;
	int *face2Sulci;

	face2Sulci = (int *)malloc(sizeof(int)*out->surface.faceNum);
	memset(face2Sulci,0,sizeof(int)*out->surface.faceNum);

	for(i=0; i<out->numberOfSulci; i++)
	{
		for(j=0; j<out->sulci[i].n; j++)
		{
			face2Sulci[out->sulci[i].faceID[j]]	++;	
		}
	}
	
	for(i=0; i<out->numberOfSulci; i++)
	{
		if(out->sulci[i].isCenter[out->sulci[i].n-1] && face2Sulci[out->sulci[i].faceID[out->sulci[i].n-1]] == 1)
		{
			out->sulci[i].n--;
			printf("remove a single center point in sulci %d\n",i);
		}
	}

	free(face2Sulci);

	return;
}

void AddNewSulciFromOldSulci(SulciTrackingOut *out, SulciTrackingOut *newOut, int sulciLabel, int start, int end)
{
	
	if(newOut->sulci)
		newOut->sulci = (Sulci *)realloc(newOut->sulci,sizeof(Sulci)*(newOut->numberOfSulci+1));
	else
		newOut->sulci = (Sulci *)malloc(sizeof(Sulci)*(newOut->numberOfSulci+1));
			
	// allocate memory
	newOut->sulci[newOut->numberOfSulci].n = end - start + 1;
	newOut->sulci[newOut->numberOfSulci].curv1 = (float *)malloc(sizeof(float)*newOut->sulci[newOut->numberOfSulci].n);
	newOut->sulci[newOut->numberOfSulci].vertex = (Fvector3d *)malloc(sizeof(Fvector3d)*newOut->sulci[newOut->numberOfSulci].n);
	newOut->sulci[newOut->numberOfSulci].pdir2 = (Fvector3d *)malloc(sizeof(Fvector3d)*newOut->sulci[newOut->numberOfSulci].n);
	newOut->sulci[newOut->numberOfSulci].segmentID = (int *)malloc(sizeof(int)*newOut->sulci[newOut->numberOfSulci].n);
	newOut->sulci[newOut->numberOfSulci].faceID = (int *)malloc(sizeof(int)*newOut->sulci[newOut->numberOfSulci].n);
	newOut->sulci[newOut->numberOfSulci].vertexID = (int *)malloc(sizeof(int)*newOut->sulci[newOut->numberOfSulci].n);
	newOut->sulci[newOut->numberOfSulci].edgeID = (int *)malloc(sizeof(int)*newOut->sulci[newOut->numberOfSulci].n);
	newOut->sulci[newOut->numberOfSulci].isJunction = (bool *)malloc(sizeof(bool)*newOut->sulci[newOut->numberOfSulci].n);
	newOut->sulci[newOut->numberOfSulci].isCenter = (bool *)malloc(sizeof(bool)*newOut->sulci[newOut->numberOfSulci].n);

	float len = 0.0;
	for(int j=0; j<newOut->sulci[newOut->numberOfSulci].n; j++)
	{
		newOut->sulci[newOut->numberOfSulci].curv1[j] = out->sulci[sulciLabel].curv1[start+j];
		
		newOut->sulci[newOut->numberOfSulci].vertex[j].x = out->sulci[sulciLabel].vertex[start+j].x;
		newOut->sulci[newOut->numberOfSulci].vertex[j].y = out->sulci[sulciLabel].vertex[start+j].y;
		newOut->sulci[newOut->numberOfSulci].vertex[j].z = out->sulci[sulciLabel].vertex[start+j].z;

		newOut->sulci[newOut->numberOfSulci].pdir2[j].x = out->sulci[sulciLabel].pdir2[start+j].x;
		newOut->sulci[newOut->numberOfSulci].pdir2[j].y = out->sulci[sulciLabel].pdir2[start+j].y;
		newOut->sulci[newOut->numberOfSulci].pdir2[j].z = out->sulci[sulciLabel].pdir2[start+j].z;

		newOut->sulci[newOut->numberOfSulci].segmentID[j] = out->sulci[sulciLabel].segmentID[start+j];
		newOut->sulci[newOut->numberOfSulci].faceID[j] = out->sulci[sulciLabel].faceID[start+j];
		newOut->sulci[newOut->numberOfSulci].edgeID[j] = out->sulci[sulciLabel].edgeID[start+j];
		newOut->sulci[newOut->numberOfSulci].vertexID[j] = out->sulci[sulciLabel].vertexID[start+j];
		newOut->sulci[newOut->numberOfSulci].isJunction[j] = out->sulci[sulciLabel].isJunction[start+j];
		newOut->sulci[newOut->numberOfSulci].isCenter[j] = out->sulci[sulciLabel].isCenter[start+j];

		if(j > 0)
			len += DistanceBetweenTwoVertices(newOut->sulci[newOut->numberOfSulci].vertex[j-1], newOut->sulci[newOut->numberOfSulci].vertex[j]);
	}
	
	newOut->sulci[newOut->numberOfSulci].color = out->sulci[sulciLabel].color;
	//newOut->sulci[newOut->numberOfSulci].junction = out->sulci[sulciLabel].junction;
	newOut->sulci[newOut->numberOfSulci].junction = 0;
	newOut->sulci[newOut->numberOfSulci].len = len;

	newOut->numberOfSulci++;
	
	return;
}

void CombineAdjacentSulci(SulciTrackingOut *out, Segment *segment, SulciTrackingOut *newOut)
{
	Surface *surface; // a pointer defined for convenience, points to the surface storage space 
	int i, j, k, m, n;
	float lenThresh;
	int *edgeSulci;
	int *edgeOrderSulci;

	int **vertexSulci;
	float **costSulci;

	surface = &(out->surface);

	// record sulci label in each edge
	edgeSulci = (int *)malloc(sizeof(int)*surface->edgeNum);
	
	// record the order in sulci in each edge
	edgeOrderSulci = (int *)malloc(sizeof(int)*surface->edgeNum);

	// record the cost of connecting two sulci
	costSulci = Falloc2d(out->numberOfSulci, out->numberOfSulci);

	// record the vertex for connecting two sulci
	vertexSulci = Ialloc2d(out->numberOfSulci, out->numberOfSulci);


	for(i=0; i<out->numberOfSulci; i++)
	for(j=0; j<out->numberOfSulci; j++)
	{
		vertexSulci[i][j] = -1;
		costSulci[i][j] = 1000.0;
	}

    for(i=0; i<surface->edgeNum; i++)
	{
		edgeSulci[i] = -1;
		edgeOrderSulci[i] = -1;
	}

	for(i=0; i<out->numberOfSulci; i++)
	{
		for(j=0; j<out->sulci[i].n; j++)
		{
			if(out->sulci[i].edgeID[j] >=0)
			{
				edgeSulci[out->sulci[i].edgeID[j]] = i;
				edgeOrderSulci[out->sulci[i].edgeID[j]] = j;
			}
		}
	}


	int *edgeSelect;
	int numSelect = 0;

	edgeSelect = (int *)malloc(sizeof(int)*surface->numMaxNeighbors);

	for(i=0; i<surface->vertexNum; i++)
	{
		numSelect = 0;
		for(j=0; j<surface->connectedges[i].size(); j++)
		{
			if(edgeSulci[surface->connectedges[i][j]] >= 0)
			{
				edgeSelect[numSelect] = surface->connectedges[i][j];
				numSelect++;
			}
		}
		
		if(numSelect >=2)
		{
			for(j=1; j<numSelect; j++)
			{
				if(out->sulci[edgeSulci[edgeSelect[0]]].color != out->sulci[edgeSulci[edgeSelect[j]]].color)
				{
					float cost = 0.0;

					for(k=0; k<numSelect; k++)
					{
						cost += DistanceBetweenTwoVertices(surface->vertex[i], out->sulci[edgeSulci[edgeSelect[k]]].vertex[edgeOrderSulci[edgeSelect[k]]]);
					}
					cost /= numSelect;

					if(cost < costSulci[edgeSulci[edgeSelect[0]]][edgeSulci[edgeSelect[j]]])
					{
						costSulci[edgeSulci[edgeSelect[0]]][edgeSulci[edgeSelect[j]]] = cost;
						costSulci[edgeSulci[edgeSelect[j]]][edgeSulci[edgeSelect[0]]] = cost;
						vertexSulci[edgeSulci[edgeSelect[0]]][edgeSulci[edgeSelect[j]]] = i;
						vertexSulci[edgeSulci[edgeSelect[j]]][edgeSulci[edgeSelect[0]]] = i;
					}
					break;
				}
			}
		}
	}

	for(m=0; m<out->numberOfSulci; m++)
	for(n=m+1; n<out->numberOfSulci; n++)
	{
		if(vertexSulci[m][n] >= 0)
		{
			i = vertexSulci[m][n];
			numSelect = 0;
			for(j=0; j<surface->connectedges[i].size(); j++)
			{
				if(edgeSulci[surface->connectedges[i][j]] >= 0)
				{
					edgeSelect[numSelect] = surface->connectedges[i][j];
					numSelect++;
				}
			}

			for(j=0; j<numSelect; j++)
			{
				out->sulci[edgeSulci[edgeSelect[j]]].junction = 1;
				out->sulci[edgeSulci[edgeSelect[j]]].isJunction[edgeOrderSulci[edgeSelect[j]]] = true;
				out->sulci[edgeSulci[edgeSelect[j]]].isCenter[edgeOrderSulci[edgeSelect[j]]] = false;
				
				out->sulci[edgeSulci[edgeSelect[j]]].vertex[edgeOrderSulci[edgeSelect[j]]].x = surface->vertex[i].x;
				out->sulci[edgeSulci[edgeSelect[j]]].vertex[edgeOrderSulci[edgeSelect[j]]].y = surface->vertex[i].y;
				out->sulci[edgeSulci[edgeSelect[j]]].vertex[edgeOrderSulci[edgeSelect[j]]].z = surface->vertex[i].z;
				
				out->sulci[edgeSulci[edgeSelect[j]]].pdir2[edgeOrderSulci[edgeSelect[j]]].x = surface->pdir2[i].x;
				out->sulci[edgeSulci[edgeSelect[j]]].pdir2[edgeOrderSulci[edgeSelect[j]]].y = surface->pdir2[i].y;
				out->sulci[edgeSulci[edgeSelect[j]]].pdir2[edgeOrderSulci[edgeSelect[j]]].z = surface->pdir2[i].z;

				out->sulci[edgeSulci[edgeSelect[j]]].edgeID[edgeOrderSulci[edgeSelect[j]]] = -1;
				out->sulci[edgeSulci[edgeSelect[j]]].faceID[edgeOrderSulci[edgeSelect[j]]] = -1;
				out->sulci[edgeSulci[edgeSelect[j]]].segmentID[edgeOrderSulci[edgeSelect[j]]] = -1;
				out->sulci[edgeSulci[edgeSelect[j]]].vertexID[edgeOrderSulci[edgeSelect[j]]] = i;				
			}
		}
	}

	int count;

	// reallocte color for each sulci
	for(i=0; i<out->numberOfSulci; i++)
	for(j=i+1; j<out->numberOfSulci; j++)
	{
		if(vertexSulci[i][j] >= 0)
		{
			if(out->sulci[i].color != out->sulci[j].color)
			{
				if(out->color2Sulci[out->sulci[i].color].size() >= out->color2Sulci[out->sulci[j].color].size())
				{
					count = out->color2Sulci[out->sulci[j].color].size();

					for(k=0; k<count; k++)
					{
						out->color2Sulci[out->sulci[i].color].push_back(out->color2Sulci[out->sulci[j].color][k]);
						out->color2Sulci[out->sulci[j].color].pop_back();
						out->sulci[out->color2Sulci[out->sulci[j].color][k]].color = out->sulci[i].color;
					}
				}
				else
				{
					count = out->color2Sulci[out->sulci[i].color].size();

					for(k=0; k<count; k++)
					{
						out->color2Sulci[out->sulci[j].color].push_back(out->color2Sulci[out->sulci[i].color][k]);
						out->color2Sulci[out->sulci[i].color].pop_back();
						out->sulci[out->color2Sulci[out->sulci[i].color][k]].color = out->sulci[j].color;					
					}
				}
			}
		}
	}

/*
	// relabel sulci color
	int *visited;
	int *indexTable;
	int neighborsFound;
	int label, index, iTmp;
	
	visited = (int *)malloc(sizeof(int)*out->numberOfSulci);
	indexTable = (int *)malloc(sizeof(int)*out->numberOfSulci);

	for(i=0; i<out->numberOfSulci; i++)
		visited[i] = -1;

	label = -1;
	for(i=0; i<out->numberOfSulci; i++)
	{
		if(out->sulci[i].junction > 0 && visited[i] < 0)
		{
			visited[i] = ++label;
			out->sulci[i].color = label;

			index = 0;
			iTmp = i;
			indexTable[index] = iTmp;

			while(1)
			{
				neighborsFound	= 0;

				for(j=0; j<out->numberOfSulci; j++)
				{
					if( vertexSulci[iTmp][j] >= 0 && out->sulci[j].junction > 0 && visited[j] < 0)
					{
						visited[j] = label;
						out->sulci[j].color = label;

						++index;
						indexTable[index] = j;
						neighborsFound = 1;
					}
				}

				if(neighborsFound)
				{
					iTmp = indexTable[index];
				}
				else if (index > 1)
				{
					--index;
					iTmp = indexTable[index];
				}
				else
					break; 
			}
		}
		else
		{
			visited[i] = ++label;
			out->sulci[n].color = label;
		}			
	}

	free(visited);
	free(indexTable);
*/

	// generate new sulci and prune small branch
	float length;
	int start, end;

	lenThresh = 8.0;

	newOut->numberOfSulci = 0;
	newOut->sulci = NULL;

	for(i=0; i<out->numberOfSulci; i++)
	{
		if(out->sulci[i].junction == 1)
		{
			start = 0;
			length = 0.0;
			
			for(j=start+1; j<out->sulci[i].n; j++)
			{
				if(out->sulci[i].isJunction[j] || j == out->sulci[i].n-1)
				{
					end = j;
					
					if(length > lenThresh)
						AddNewSulciFromOldSulci(out, newOut, i, start, end);
	
					start = end;
					length = 0;
				}
				length += DistanceBetweenTwoVertices(out->sulci[i].vertex[j-1], out->sulci[i].vertex[j]);
			}
		}
		else
			AddNewSulciFromOldSulci(out, newOut, i, 0, out->sulci[i].n-1);
	}

	free(edgeSelect);
	free(edgeSulci);
	free(edgeOrderSulci);
	Ifree2d(vertexSulci,out->numberOfSulci);
	Ffree2d(costSulci,out->numberOfSulci);

	return;
}

void FindNearstVertex(Surface *surface, int edgeID, Fvector3d halfVertex0, Fvector3d halfVertex1, int *vertexID)
{
	Fvector3d tempVector0, tempVector1, tempVector2;
	Fvector3d *vertex1, *vertex2;

	tempVector0 = Fvector3dMINUSFvector3d( halfVertex1,  halfVertex0);

	vertex1 = &(surface->vertex[surface->edges[edgeID].vertexID[0]]);

	vertex2 = &(surface->vertex[surface->edges[edgeID].vertexID[1]]);

	tempVector1 = Fvector3dMINUSFvector3d( *vertex1,  halfVertex0);
	
	tempVector2 = Fvector3dMINUSFvector3d( *vertex2,  halfVertex0);

	if( Fvector3dDOTFvector3d(tempVector0, tempVector1) < 0)
		*vertexID = surface->edges[edgeID].vertexID[0];
	else
		*vertexID = surface->edges[edgeID].vertexID[1];
	
	return;
}



void AddSulciSegment(Surface *surface, SulciTrackingOut *out, Segment *segment, 
					 int *faceSegment, int *faceSulci, int sulcusLabel, 
					 Fvector3d startVertex, int startVertexID, int endVertexID, VtxOnFace *pathVtxList, int pathLen)
{

	Fvector3d V1, V2, V;
	double l;

	if(out->sulci)
		out->sulci = (Sulci*)realloc(out->sulci,sizeof(Sulci)*(out->numberOfSulci+1));
	else
		out->sulci = (Sulci*)malloc(sizeof(Sulci)*(out->numberOfSulci+1));
			
	out->sulci[out->numberOfSulci].n = pathLen + 2;
	out->sulci[out->numberOfSulci].vertex = (Fvector3d *)malloc(sizeof(Fvector3d)*out->sulci[out->numberOfSulci].n);
	out->sulci[out->numberOfSulci].curv1 = (float *)malloc(sizeof(float)*out->sulci[out->numberOfSulci].n);

	out->sulci[out->numberOfSulci].curv1[0] = 0.0;
/*
	out->sulci[out->numberOfSulci].vertex[0].x = surface->vertex[startVertexID].x;
	out->sulci[out->numberOfSulci].vertex[0].y = surface->vertex[startVertexID].y;
	out->sulci[out->numberOfSulci].vertex[0].z = surface->vertex[startVertexID].z;
*/
	out->sulci[out->numberOfSulci].vertex[0].x = startVertex.x;
	out->sulci[out->numberOfSulci].vertex[0].y = startVertex.y;
	out->sulci[out->numberOfSulci].vertex[0].z = startVertex.z;

	for(int i = 0; i< pathLen; i++)
	{
		if(pathVtxList[i].lambda==0.0)
		{
			V.x = surface->vertex[pathVtxList[i].ID].x;
			V.y = surface->vertex[pathVtxList[i].ID].y;
			V.z = surface->vertex[pathVtxList[i].ID].z;
		}
		else
		{
			l = pathVtxList[i].lambda;
			
			V1.x = surface->vertex[surface->edges[pathVtxList[i].ID].vertexID[0]].x;
			V1.y = surface->vertex[surface->edges[pathVtxList[i].ID].vertexID[0]].y;
			V1.z = surface->vertex[surface->edges[pathVtxList[i].ID].vertexID[0]].z;
			V2.x = surface->vertex[surface->edges[pathVtxList[i].ID].vertexID[1]].x;
			V2.y = surface->vertex[surface->edges[pathVtxList[i].ID].vertexID[1]].y;
			V2.z = surface->vertex[surface->edges[pathVtxList[i].ID].vertexID[1]].z;
			
			V.x = V1.x+l*(V2.x-V1.x);
			V.y = V1.y+l*(V2.y-V1.y);
			V.z = V1.z+l*(V2.z-V1.z);
		 }

		out->sulci[out->numberOfSulci].curv1[pathLen-i] = 0.0;
		out->sulci[out->numberOfSulci].vertex[pathLen-i].x = V.x;
		out->sulci[out->numberOfSulci].vertex[pathLen-i].y = V.y;
		out->sulci[out->numberOfSulci].vertex[pathLen-i].z = V.z;

	}
			
	float distS, distE, minDist;
	Fvector3d *selectVertex;

	minDist = 1000.0;
	for(int i=0; i<surface->adjacentfaces[endVertexID].size(); i++)
	{
		if(faceSulci[surface->adjacentfaces[endVertexID][i]] >= 0)// && faceSulci[surface->adjacentfaces[endVertexID][i]] != sulcusLabel)
		{
			distS = DistanceBetweenTwoVertices(surface->vertex[endVertexID], segment->start[faceSegment[surface->adjacentfaces[endVertexID][i]]]);
			distE = DistanceBetweenTwoVertices(surface->vertex[endVertexID], segment->end[faceSegment[surface->adjacentfaces[endVertexID][i]]]);

			if(distS < minDist)
			{
				minDist = distS;
				selectVertex = &(segment->start[faceSegment[surface->adjacentfaces[endVertexID][i]]]);
			}
			if(distE < minDist)
			{
				minDist = distE;
				selectVertex = &(segment->end[faceSegment[surface->adjacentfaces[endVertexID][i]]]);
			}
		}
	}

	if(selectVertex->x >= 0 && selectVertex->y >= 0 && selectVertex->z >= 0)
	{
		out->sulci[out->numberOfSulci].vertex[pathLen+1].x = selectVertex->x;
		out->sulci[out->numberOfSulci].vertex[pathLen+1].y = selectVertex->y;
		out->sulci[out->numberOfSulci].vertex[pathLen+1].z = selectVertex->z;
	}
/*
	out->sulci[out->numberOfSulci].vertex[pathLen+1].x = surface->vertex[endVertexID].x;
	out->sulci[out->numberOfSulci].vertex[pathLen+1].y = surface->vertex[endVertexID].y;
	out->sulci[out->numberOfSulci].vertex[pathLen+1].z = surface->vertex[endVertexID].z;
*/
	out->sulci[out->numberOfSulci].curv1[pathLen+1] = 0.0;

	out->numberOfSulci++;

	return;
}


void ConnectInterruptedSulci(Surface *surface, Segment *segment, SulciTrackingOut *out, SulciTrackingOut *outTest)
{  
	int *faceSulci;
	int *faceSegment;
	double *marchingSpeed;
	int *VtxID;
	VtxOnFace *pathVtxList;
	int pathLen;
	int N;

	int i, j;
	int faceID;
	int startVertexID, endVertexID;
	float geodesicThresh;

	// record sulci label in each face
	faceSulci = (int *)malloc(sizeof(int)*surface->faceNum);

	// record segment ID in each face
	faceSegment = (int *)malloc(sizeof(int)*surface->faceNum);

	marchingSpeed = (double *)malloc(sizeof(double)*surface->vertexNum);

	pathVtxList = (VtxOnFace *)malloc(sizeof(VtxOnFace)*1000);

	VtxID = (int *)malloc(sizeof(int)*100);

	for(i=0; i<surface->faceNum; i++)
	{
		faceSulci[i] = -1;
		faceSegment[i] = -1;
	}

	// uniform marching speed
	for(i=0; i<surface->vertexNum; i++)
		marchingSpeed[i] = 1.0;

	for(i=0; i<out->numberOfSulci; i++)
	{
		for(j=0; j<out->sulci[i].n; j++)
		{
			faceSulci[out->sulci[i].faceID[j]] = i;
			faceSegment[out->sulci[i].faceID[j]] = out->sulci[i].segmentID[j];
		}
	}

	geodesicThresh = 8.0;

	outTest->numberOfSulci = 0;
	outTest->sulci = (Sulci*)malloc(sizeof(Sulci)*(outTest->numberOfSulci+1));

	for(i=0; i<out->numberOfSulci; i++)
	{

		startVertexID = -1;

		FindNearstVertex(surface, out->sulci[i].edgeID[0], out->sulci[i].vertex[0], out->sulci[i].vertex[1], &startVertexID);

		if(startVertexID >= 0 && startVertexID < surface->vertexNum)
		{
			endVertexID = -1;
			SearchConnectingVertexUsingFastMarching(surface, startVertexID, marchingSpeed, faceSulci, geodesicThresh, i, &endVertexID);

			if( endVertexID >=0 && endVertexID < surface->vertexNum)
			{
				N = 2;
				VtxID[0] = startVertexID;
				VtxID[1] = endVertexID;
			
				PathFinder(surface, marchingSpeed, VtxID, N, pathVtxList, &pathLen);
				//pathLen = 0;
				AddSulciSegment(surface, outTest, segment, faceSegment, faceSulci, i, out->sulci[i].vertex[0], startVertexID, endVertexID, pathVtxList, pathLen);
			}
		}

		
		startVertexID = -1;

		FindNearstVertex(surface, out->sulci[i].edgeID[out->sulci[i].n-1], out->sulci[i].vertex[out->sulci[i].n-1], out->sulci[i].vertex[out->sulci[i].n-2], &startVertexID);

		if(startVertexID >= 0 && startVertexID < surface->vertexNum)
		{
			endVertexID = -1;
			SearchConnectingVertexUsingFastMarching(surface, startVertexID, marchingSpeed, faceSulci, geodesicThresh, i, &endVertexID);
			
			if( endVertexID >=0 && endVertexID < surface->vertexNum)
			{
				N = 2;
				VtxID[0] = startVertexID;
				VtxID[1] = endVertexID;
			
				PathFinder(surface, marchingSpeed, VtxID, N, pathVtxList, &pathLen);
				//pathLen = 0;
				AddSulciSegment(surface, outTest, segment, faceSegment, faceSulci, i, out->sulci[i].vertex[out->sulci[i].n-1], startVertexID, endVertexID, pathVtxList, pathLen);
			}
		}
		printf("i = %d\n",i);
	}

	free(faceSulci);
	free(faceSegment);
	free(marchingSpeed);

	return;
}


void FindStart(Surface *surface, int edgeID, Fvector3d halfVertex, VtxOnFace *start)
{
	float dist1, dist2;

	dist1 = DistanceBetweenTwoVertices(surface->vertex[surface->edges[edgeID].vertexID[0]], surface->vertex[surface->edges[edgeID].vertexID[1]]);
	dist2 = DistanceBetweenTwoVertices(halfVertex, surface->vertex[surface->edges[edgeID].vertexID[0]]);
	
	start->lambda = dist2/dist1;
	if(start->lambda < 0.000001)
	{
		start->lambda = 0.0;
		start->ID = surface->edges[edgeID].vertexID[0];
	}
	else
        start->ID = edgeID;


	return;
}

Fvector3d InterpolateHalfVertex(Fvector3d V1, Fvector3d V2, float lambda)
{
	Fvector3d V;

	V.x = V1.x +lambda*(V2.x-V1.x);
	V.y = V1.y +lambda*(V2.y-V1.y);
	V.z = V1.z +lambda*(V2.z-V1.z);

	return V;
}

void AddSulciSegmentHalfVertex(Surface *surface, SulciTrackingOut *out, VtxOnFace start, VtxOnFace end, VtxOnFace *pathVtxList, int pathLen)
{

	Fvector3d V1, V2, V;
	double l;

	if(out->sulci)
		out->sulci = (Sulci*)realloc(out->sulci,sizeof(Sulci)*(out->numberOfSulci+1));
	else
		out->sulci = (Sulci*)malloc(sizeof(Sulci)*(out->numberOfSulci+1));
			
	out->sulci[out->numberOfSulci].n = pathLen + 2;
	out->sulci[out->numberOfSulci].vertex = (Fvector3d *)malloc(sizeof(Fvector3d)*out->sulci[out->numberOfSulci].n);
	out->sulci[out->numberOfSulci].curv1 = (float *)malloc(sizeof(float)*out->sulci[out->numberOfSulci].n);

	out->sulci[out->numberOfSulci].vertex[0] = InterpolateHalfVertex(surface->vertex[surface->edges[start.ID].vertexID[0]], surface->vertex[surface->edges[start.ID].vertexID[1]], start.lambda);

	out->sulci[out->numberOfSulci].curv1[0] = 0.0;

	for(int i = 0; i< pathLen; i++)
	{
		if(pathVtxList[i].lambda==0.0)
		{
			V.x = surface->vertex[pathVtxList[i].ID].x;
			V.y = surface->vertex[pathVtxList[i].ID].y;
			V.z = surface->vertex[pathVtxList[i].ID].z;
		}
		else
		{
			V = InterpolateHalfVertex(surface->vertex[surface->edges[pathVtxList[i].ID].vertexID[0]], surface->vertex[surface->edges[pathVtxList[i].ID].vertexID[1]], pathVtxList[i].lambda);
		 }

		out->sulci[out->numberOfSulci].curv1[pathLen-i] = 0.0;
		out->sulci[out->numberOfSulci].vertex[pathLen-i].x = V.x;
		out->sulci[out->numberOfSulci].vertex[pathLen-i].y = V.y;
		out->sulci[out->numberOfSulci].vertex[pathLen-i].z = V.z;
	}
			
	out->sulci[out->numberOfSulci].vertex[pathLen+1] = InterpolateHalfVertex(surface->vertex[surface->edges[end.ID].vertexID[0]], surface->vertex[surface->edges[end.ID].vertexID[1]], end.lambda);

	out->sulci[out->numberOfSulci].curv1[pathLen+1] = 0.0;

	out->numberOfSulci++;

	return;
}

void AddTwoSulciSegmentHalfVertex(Surface *surface, SulciTrackingOut *out, SulciTrackingOut *newOut, int sulcusLabel, VtxOnFace end1, VtxOnFace *pathVtxList1, int pathLen1, VtxOnFace end2, VtxOnFace *pathVtxList2, int pathLen2)
{
	int i;
	Fvector3d V;

	if(pathLen1 == 0 && pathLen2 == 0)
	{
		newOut->sulci[sulcusLabel].n = out->sulci[sulcusLabel].n;
		newOut->sulci[sulcusLabel].color = out->sulci[sulcusLabel].color;
		newOut->sulci[sulcusLabel].curv1 = (float *)malloc(sizeof(float)*out->sulci[sulcusLabel].n);
		newOut->sulci[sulcusLabel].vertex = (Fvector3d *)malloc(sizeof(Fvector3d)*out->sulci[sulcusLabel].n);
		
		newOut->sulci[sulcusLabel].vertexID = (int *)malloc(sizeof(int)*out->sulci[sulcusLabel].n);
		newOut->sulci[sulcusLabel].edgeID = (int *)malloc(sizeof(int)*out->sulci[sulcusLabel].n);
		newOut->sulci[sulcusLabel].faceID = (int *)malloc(sizeof(int)*out->sulci[sulcusLabel].n);

		for(i=0; i<out->sulci[sulcusLabel].n; i++)
		{
			newOut->sulci[sulcusLabel].curv1[i] = out->sulci[sulcusLabel].curv1[i];

			newOut->sulci[sulcusLabel].vertex[i].x = out->sulci[sulcusLabel].vertex[i].x;
			newOut->sulci[sulcusLabel].vertex[i].y = out->sulci[sulcusLabel].vertex[i].y;
			newOut->sulci[sulcusLabel].vertex[i].z = out->sulci[sulcusLabel].vertex[i].z;

			newOut->sulci[sulcusLabel].vertexID[i] = out->sulci[sulcusLabel].vertexID[i];
			newOut->sulci[sulcusLabel].edgeID[i] = out->sulci[sulcusLabel].edgeID[i];
			newOut->sulci[sulcusLabel].faceID[i] = out->sulci[sulcusLabel].faceID[i];
		}
	}
	else
	{
		newOut->sulci[sulcusLabel].n = out->sulci[sulcusLabel].n + pathLen1 + pathLen2 + (pathLen1 > 0) + (pathLen2 > 0);
		newOut->sulci[sulcusLabel].color = out->sulci[sulcusLabel].color;
		newOut->sulci[sulcusLabel].curv1 = (float *)malloc(sizeof(float)*newOut->sulci[sulcusLabel].n);
		newOut->sulci[sulcusLabel].vertex = (Fvector3d *)malloc(sizeof(Fvector3d)*newOut->sulci[sulcusLabel].n);
		
		newOut->sulci[sulcusLabel].vertexID = (int *)malloc(sizeof(int)*newOut->sulci[sulcusLabel].n);
		newOut->sulci[sulcusLabel].edgeID = (int *)malloc(sizeof(int)*newOut->sulci[sulcusLabel].n);
		newOut->sulci[sulcusLabel].faceID = (int *)malloc(sizeof(int)*newOut->sulci[sulcusLabel].n);
		
		if(pathLen1 > 0)
		{
			if(end1.lambda > 0.0 )
			{
				V = InterpolateHalfVertex(surface->vertex[surface->edges[end1.ID].vertexID[0]], 
					surface->vertex[surface->edges[end1.ID].vertexID[1]], end1.lambda);

				newOut->sulci[sulcusLabel].vertexID[0] = -1;
				newOut->sulci[sulcusLabel].edgeID[0] = end1.ID;
				newOut->sulci[sulcusLabel].faceID[0] = -1;
			}
			else if(end1.lambda == 0.0)
			{	
				V.x = surface->vertex[end1.ID].x;
				V.y = surface->vertex[end1.ID].y;
				V.z = surface->vertex[end1.ID].z;
			
				newOut->sulci[sulcusLabel].vertexID[0] = end1.ID;
				newOut->sulci[sulcusLabel].edgeID[0] = -1;
				newOut->sulci[sulcusLabel].faceID[0] = -1;
			}

			newOut->sulci[sulcusLabel].curv1[0] =  0.0;
			newOut->sulci[sulcusLabel].vertex[0].x = V.x;
			newOut->sulci[sulcusLabel].vertex[0].y = V.y;
			newOut->sulci[sulcusLabel].vertex[0].z = V.z;
			

			for(i=0; i<pathLen1; i++)
			{
				if(pathVtxList1[i].lambda==0.0)
				{
					V.x = surface->vertex[pathVtxList1[i].ID].x;
					V.y = surface->vertex[pathVtxList1[i].ID].y;
					V.z = surface->vertex[pathVtxList1[i].ID].z;

					newOut->sulci[sulcusLabel].vertexID[i+1] = pathVtxList1[i].ID;
					newOut->sulci[sulcusLabel].edgeID[i+1] = -1;
					newOut->sulci[sulcusLabel].faceID[i+1] = -1;
				}
				else
				{
					V = InterpolateHalfVertex(surface->vertex[surface->edges[pathVtxList1[i].ID].vertexID[0]], 
						surface->vertex[surface->edges[pathVtxList1[i].ID].vertexID[1]], pathVtxList1[i].lambda);

					newOut->sulci[sulcusLabel].vertexID[i+1] = -1;
					newOut->sulci[sulcusLabel].edgeID[i+1] = pathVtxList1[i].ID;
					newOut->sulci[sulcusLabel].faceID[i+1] = -1;
				}
	
				newOut->sulci[sulcusLabel].curv1[i+1] =  0.0;
				newOut->sulci[sulcusLabel].vertex[i+1].x = V.x;
				newOut->sulci[sulcusLabel].vertex[i+1].y = V.y;
				newOut->sulci[sulcusLabel].vertex[i+1].z = V.z;
			}

			for(i=0; i<out->sulci[sulcusLabel].n; i++)
			{
				newOut->sulci[sulcusLabel].curv1[i+pathLen1+1] = out->sulci[sulcusLabel].curv1[i];
				
				newOut->sulci[sulcusLabel].vertex[i+pathLen1+1].x = out->sulci[sulcusLabel].vertex[i].x;
				newOut->sulci[sulcusLabel].vertex[i+pathLen1+1].y = out->sulci[sulcusLabel].vertex[i].y;
				newOut->sulci[sulcusLabel].vertex[i+pathLen1+1].z = out->sulci[sulcusLabel].vertex[i].z;

				newOut->sulci[sulcusLabel].vertexID[i+pathLen1+1] = out->sulci[sulcusLabel].vertexID[i];
				newOut->sulci[sulcusLabel].edgeID[i+pathLen1+1] = out->sulci[sulcusLabel].edgeID[i];
				newOut->sulci[sulcusLabel].faceID[i+pathLen1+1] = out->sulci[sulcusLabel].faceID[i];
			}
		}
		else
		{
			for(i=0; i<out->sulci[sulcusLabel].n; i++)
			{
				newOut->sulci[sulcusLabel].curv1[i] = out->sulci[sulcusLabel].curv1[i];
				newOut->sulci[sulcusLabel].vertex[i].x = out->sulci[sulcusLabel].vertex[i].x;
				newOut->sulci[sulcusLabel].vertex[i].y = out->sulci[sulcusLabel].vertex[i].y;
				newOut->sulci[sulcusLabel].vertex[i].z = out->sulci[sulcusLabel].vertex[i].z;

				newOut->sulci[sulcusLabel].vertexID[i] = out->sulci[sulcusLabel].vertexID[i];
				newOut->sulci[sulcusLabel].edgeID[i] = out->sulci[sulcusLabel].edgeID[i];
				newOut->sulci[sulcusLabel].faceID[i] = out->sulci[sulcusLabel].faceID[i];
			}
		}

		if(pathLen2 > 0)
		{
			for(i=0; i<pathLen2; i++)
			{
				if(pathVtxList2[i].lambda==0.0)
				{
					V.x = surface->vertex[pathVtxList2[i].ID].x;
					V.y = surface->vertex[pathVtxList2[i].ID].y;
					V.z = surface->vertex[pathVtxList2[i].ID].z;

					newOut->sulci[sulcusLabel].vertexID[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2-i-1] = pathVtxList2[i].ID;
					newOut->sulci[sulcusLabel].edgeID[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2-i-1] = -1;
					newOut->sulci[sulcusLabel].faceID[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2-i-1] = -1;
				}
				else
				{
					V = InterpolateHalfVertex(surface->vertex[surface->edges[pathVtxList2[i].ID].vertexID[0]], 
						surface->vertex[surface->edges[pathVtxList2[i].ID].vertexID[1]], pathVtxList2[i].lambda);

					newOut->sulci[sulcusLabel].vertexID[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2-i-1] = -1;
					newOut->sulci[sulcusLabel].edgeID[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2-i-1] = pathVtxList2[i].ID;
					newOut->sulci[sulcusLabel].faceID[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2-i-1] = -1;
				}
				newOut->sulci[sulcusLabel].curv1[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2-i-1] =  0.0;
				newOut->sulci[sulcusLabel].vertex[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2-i-1].x = V.x;
				newOut->sulci[sulcusLabel].vertex[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2-i-1].y = V.y;
				newOut->sulci[sulcusLabel].vertex[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2-i-1].z = V.z;
			}
			
			if(end2.lambda > 0.0)
			{
				V = InterpolateHalfVertex(surface->vertex[surface->edges[end2.ID].vertexID[0]], 
				surface->vertex[surface->edges[end2.ID].vertexID[1]], end2.lambda);

				newOut->sulci[sulcusLabel].vertexID[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2] = -1;
				newOut->sulci[sulcusLabel].edgeID[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2] = end2.ID;
				newOut->sulci[sulcusLabel].faceID[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2] = -1;
			}
			else if(end2.lambda == 0.0)
			{
				V.x = surface->vertex[end2.ID].x;
				V.y = surface->vertex[end2.ID].y;
				V.z = surface->vertex[end2.ID].z;

				newOut->sulci[sulcusLabel].vertexID[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2] = end2.ID;
				newOut->sulci[sulcusLabel].edgeID[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2] = -1;
				newOut->sulci[sulcusLabel].faceID[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2] = -1;
			}
		
			newOut->sulci[sulcusLabel].curv1[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2] =  0.0;
			newOut->sulci[sulcusLabel].vertex[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2].x = V.x;
			newOut->sulci[sulcusLabel].vertex[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2].y = V.y;
			newOut->sulci[sulcusLabel].vertex[pathLen1+(pathLen1>0)+out->sulci[sulcusLabel].n+pathLen2].z = V.z;
		}
	}

	return;
}

void AddOneSulciToOldSulci(SulciTrackingOut *out, int sulciLabel, int start, int end)
{
	
	out->sulci = (Sulci *)realloc(out->sulci,sizeof(Sulci)*(out->numberOfSulci+1));
			
	// allocate memory
	out->sulci[out->numberOfSulci].n = end - start + 1;
	out->sulci[out->numberOfSulci].curv1 = (float *)malloc(sizeof(float)*out->sulci[out->numberOfSulci].n);
	out->sulci[out->numberOfSulci].vertex = (Fvector3d *)malloc(sizeof(Fvector3d)*out->sulci[out->numberOfSulci].n);
	out->sulci[out->numberOfSulci].faceID = (int *)malloc(sizeof(int)*out->sulci[out->numberOfSulci].n);
	out->sulci[out->numberOfSulci].vertexID = (int *)malloc(sizeof(int)*out->sulci[out->numberOfSulci].n);
	out->sulci[out->numberOfSulci].edgeID = (int *)malloc(sizeof(int)*out->sulci[out->numberOfSulci].n);
//	out->sulci[out->numberOfSulci].isJunction = (bool *)malloc(sizeof(bool)*out->sulci[out->numberOfSulci].n);

	float len = 0.0;
	for(int j=0; j<out->sulci[out->numberOfSulci].n; j++)
	{
		out->sulci[out->numberOfSulci].curv1[j] = out->sulci[sulciLabel].curv1[start+j];
		
		out->sulci[out->numberOfSulci].vertex[j].x = out->sulci[sulciLabel].vertex[start+j].x;
		out->sulci[out->numberOfSulci].vertex[j].y = out->sulci[sulciLabel].vertex[start+j].y;
		out->sulci[out->numberOfSulci].vertex[j].z = out->sulci[sulciLabel].vertex[start+j].z;

		out->sulci[out->numberOfSulci].faceID[j] = out->sulci[sulciLabel].faceID[start+j];
		out->sulci[out->numberOfSulci].edgeID[j] = out->sulci[sulciLabel].edgeID[start+j];
		out->sulci[out->numberOfSulci].vertexID[j] = out->sulci[sulciLabel].vertexID[start+j];
//		out->sulci[out->numberOfSulci].isJunction[j] = out->sulci[sulciLabel].isJunction[start+j];

		if(j > 0)
			len += DistanceBetweenTwoVertices(out->sulci[out->numberOfSulci].vertex[j-1], out->sulci[out->numberOfSulci].vertex[j]);
	}
	
	out->sulci[out->numberOfSulci].color = out->sulci[sulciLabel].color;
	out->sulci[out->numberOfSulci].junction = 0;
	out->sulci[out->numberOfSulci].len = len;

	out->numberOfSulci++;
	
	return;
}


void InterruptSulci(SulciTrackingOut *out, SulciTrackingOut *newOut)
{
	int i, j;

	int junction;
	
	vector< vector<int> > junctionPoint; 	

	junctionPoint.resize(out->numberOfSulci);

	for(i=0; i<out->numberOfSulci; i++)
	{
		if(out->sulci[i].junction)
		{
			for(j=0; j<out->sulci[i].n; j++)
			{
				if(out->sulci[i].isJunction[j])
				{
					if(out->sulci[i].edgeID[j] >=0 )
						junction = out->sulci[i].edgeID[j];
					else if(out->sulci[i].vertexID[j] >=0 )
						junction = out->sulci[i].vertexID[j];
					else
						junction = out->sulci[i].faceID[j];
					
					junctionPoint[i].push_back(junction);
				}
			}
		}
	}
	
	for(i=0; i<out->numberOfSulci; i++)
		printf("sulci %d: junction %d\n",i,junctionPoint[i].size());

	int k, n;
	int start, end;
	bool findJunction;

	for(i=0; i<out->numberOfSulci; i++)
	{
		if(out->sulci[i].junction)
		{
			start  = 0;
			for(k=start+1; k<newOut->sulci[i].n; k++)
			{
				findJunction = false;
				for(n=0; n<junctionPoint[i].size(); n++)
				{
					if( (newOut->sulci[i].edgeID[k] == junctionPoint[i][n] || newOut->sulci[i].vertexID[k] == junctionPoint[i][n]) || k == newOut->sulci[i].n-1)
					{
						end = 	k;
						findJunction = true;
						break;
					}
				}
				if(findJunction)
				{
					AddOneSulciToOldSulci(newOut, i, start, end);
					start = end;
				}				
			}
			newOut->sulci[i].n = 0;
		}
	}

	return;
}

float InterpolateCurvature(Surface *surface, VtxOnFace edge)
{
	float curvature;

	if(edge.lambda == 0.0)
		curvature = surface->curv1[edge.ID];
	else
	{
		/*
		if(surface->curv1[surface->edges[edge.ID].vertexID[0]]  > surface->curv1[surface->edges[edge.ID].vertexID[1]])
			curvature = surface->curv1[surface->edges[edge.ID].vertexID[0]];
		else
			curvature = surface->curv1[surface->edges[edge.ID].vertexID[1]];
		*/
		curvature = (1.0-edge.lambda)*surface->curv1[surface->edges[edge.ID].vertexID[0]] + edge.lambda*surface->curv1[surface->edges[edge.ID].vertexID[1]];
	}

	return curvature;
}

void ComputeNonUniformlMarchingSpeed(Surface *pmesh, double* marchingSpeed)
{
	float alpha = -5.0; 
	float thresh = -0.5;

	for(int i=0; i<pmesh->vertexNum; i++)
	{
		if(pmesh->curv1[i] < thresh)
			marchingSpeed[i] = 1.0;
		else if (pmesh->curv1[i] > -1.0*thresh)
			marchingSpeed[i] = exp(alpha);
		else
			marchingSpeed[i] = exp(alpha*fabs(pmesh->curv1[i] - thresh));

		if(pmesh->neighbors[i].size() != pmesh->adjacentfaces[i].size())
			marchingSpeed[i] = 0.00000001;
	}

	return;
}

void ComputeNonUniformlMarchingSpeedNew(Surface *pmesh, Segment *segment, double* marchingSpeed)
{
	int i;
	float alpha = -5.0; 
	float thresh = -0.5;
	float belta = 0.8;

	char *vertexType;

	vertexType = (char *)malloc(sizeof(char)*pmesh->vertexNum);

	memset(vertexType,0,sizeof(char)*pmesh->vertexNum);

	for(i=0; i<segment->n; i++)
	{
		vertexType[pmesh->faces[segment->face[i]].v[0]] = 1;
		vertexType[pmesh->faces[segment->face[i]].v[1]] = 1;
		vertexType[pmesh->faces[segment->face[i]].v[2]] = 1;
	}
	
	for(int i=0; i<pmesh->vertexNum; i++)
	{
		if(pmesh->curv1[i] < thresh)
			marchingSpeed[i] = 1.0;
		else if (pmesh->curv1[i] > -1.0*thresh)
			marchingSpeed[i] = exp(alpha);
		else
			marchingSpeed[i] = exp(alpha*fabs(pmesh->curv1[i] - thresh));

		marchingSpeed[i] = belta*marchingSpeed[i] + (1.0-belta)*vertexType[i];

		if(pmesh->neighbors[i].size() != pmesh->adjacentfaces[i].size())
			marchingSpeed[i] = 0.00000001;
	}

	free(vertexType);

	return;
}

float ComputeRatio(Fvector3d start, Fvector3d end, Fvector3d mid)
{
	return(DistanceBetweenTwoVertices(start, mid)/DistanceBetweenTwoVertices(start, end));
}

void ComputeNonUniformlMarchingSpeedNew2(Surface *pmesh, Segment *segment, double* marchingSpeed)
{
	int i;
	float alpha = -5.0; 
	float thresh = -0.5;
	float belta = 0.5;

	float *vertexValue;

	vertexValue = (float *)malloc(sizeof(float)*pmesh->vertexNum);
	
	for(i=0; i<pmesh->vertexNum; i++)
		vertexValue[i] = 0.0000001;

	for(i=0; i<segment->n; i++)
	{
		float ratio;

		ratio = ComputeRatio(pmesh->vertex[pmesh->edges[segment->edgeS[i]].vertexID[0]], pmesh->vertex[pmesh->edges[segment->edgeS[i]].vertexID[1]], segment->start[i]);

		if(vertexValue[pmesh->edges[segment->edgeS[i]].vertexID[0]] < 1-ratio )
			vertexValue[pmesh->edges[segment->edgeS[i]].vertexID[0]] = 1-ratio;
	}

	for(int i=0; i<pmesh->vertexNum; i++)
	{
		if(pmesh->curv1[i] < thresh)
			marchingSpeed[i] = 1.0;
		else if (pmesh->curv1[i] > -1.0*thresh)
			marchingSpeed[i] = exp(alpha);
		else
			marchingSpeed[i] = exp(alpha*fabs(pmesh->curv1[i] - thresh));

		marchingSpeed[i] = belta*marchingSpeed[i] + (1.0-belta)*vertexValue[i];

		if(pmesh->neighbors[i].size() != pmesh->adjacentfaces[i].size())
			marchingSpeed[i] = 0.00000001;
	}

	free(vertexValue);

	return;
}



void ConnectInterruptedSulciHalfVertex(Surface *surface, Segment *segment, SulciTrackingOut *out, SulciTrackingOut *newOut)
{  
	int *edgeSulci;
	unsigned char *edgeSulciEnd;
	Fvector3d *edgeHalfVertex;
	bool *sulciEndVisited;
	bool *sulciStartVisited;
	double *marchingSpeed;
	double *marchingSpeedNU;
	VtxOnFace *pathVtxList1, *pathVtxList2;
	int pathLen1, pathLen2;
	int N;

	int i, j;
	int faceID;
	int startVertexID, endVertexID;
	float geodesicThresh;
	float curvThresh;

	// record sulci label in each edge
	edgeSulci = (int *)malloc(sizeof(int)*surface->edgeNum);

	// record whether is sulci end in each edge
	edgeSulciEnd = (unsigned char*)malloc(sizeof(unsigned char)*surface->edgeNum);

	// record sulci position in each edge
	edgeHalfVertex = (Fvector3d *)malloc(sizeof(Fvector3d)*surface->edgeNum);

	sulciEndVisited = (bool *)malloc(sizeof(bool)*out->numberOfSulci);
	sulciStartVisited = (bool *)malloc(sizeof(bool)*out->numberOfSulci);

	marchingSpeed = (double *)malloc(sizeof(double)*surface->vertexNum);
	marchingSpeedNU = (double *)malloc(sizeof(double)*surface->vertexNum);

	pathVtxList1 = (VtxOnFace *)malloc(sizeof(VtxOnFace)*1000);
	pathVtxList2 = (VtxOnFace *)malloc(sizeof(VtxOnFace)*1000);

	for(i=0; i<surface->edgeNum; i++)
	{
		edgeSulci[i] = -1;
		edgeSulciEnd[i] = 0;
		edgeHalfVertex[i].x = -1.0;
		edgeHalfVertex[i].y = -1.0;
		edgeHalfVertex[i].z = -1.0;
	}

	for(i=0; i<out->numberOfSulci; i++)
	{
		sulciStartVisited[i] = false;
		sulciEndVisited[i] = false;
	}


	// uniform marching speed
	for(i=0; i<surface->vertexNum; i++)
	{
		marchingSpeed[i] = 1.0;

		if(surface->neighbors[i].size() != surface->adjacentfaces[i].size())
			marchingSpeed[i] = 0.00000001;
	}

	// nonuniform marching speed
	ComputeNonUniformlMarchingSpeed(surface, marchingSpeedNU);

	for(i=0; i<out->numberOfSulci; i++)
	{
		if(out->sulci[i].edgeID[0] < 0)
		{
			for(j=0; j<out->sulci[i].n; j++)
			{
				out->sulci[i].curv1[j] = out->sulci[i].curv1[j+1];
				out->sulci[i].edgeID[j] = out->sulci[i].edgeID[j+1];
				out->sulci[i].faceID[j] = out->sulci[i].faceID[j+1];
				out->sulci[i].segmentID[j] = out->sulci[i].segmentID[j+1];
				
				out->sulci[i].vertex[j].x = out->sulci[i].vertex[j+1].x;
				out->sulci[i].vertex[j].y = out->sulci[i].vertex[j+1].y;
				out->sulci[i].vertex[j].z = out->sulci[i].vertex[j+1].z;

				out->sulci[i].pdir2[j].x = out->sulci[i].pdir2[j+1].x;
				out->sulci[i].pdir2[j].y = out->sulci[i].pdir2[j+1].y;
				out->sulci[i].pdir2[j].z = out->sulci[i].pdir2[j+1].z;
			}
			out->sulci[i].n = out->sulci[i].n - 1; 
		}
		if(out->sulci[i].edgeID[out->sulci[i].n-1] < 0)
		{
			out->sulci[i].n = out->sulci[i].n - 1; 
		}
	}


	for(i=0; i<out->numberOfSulci; i++)
	{
		for(j=0; j<out->sulci[i].n; j++)
		{
			if(out->sulci[i].edgeID[j] <0 || out->sulci[i].edgeID[j] >= surface->edgeNum)
				printf("Error: i = %d, j = %d, edgeID = %d, len = %d \n", i, j, out->sulci[i].edgeID[j], out->sulci[i].n);

			// start point
			if(j==0)
				edgeSulciEnd[out->sulci[i].edgeID[j]] = 1;
			//end point
			else if(j == out->sulci[i].n-1)
				edgeSulciEnd[out->sulci[i].edgeID[j]] = 2;

			edgeSulci[out->sulci[i].edgeID[j]] = i;
			edgeHalfVertex[out->sulci[i].edgeID[j]].x = out->sulci[i].vertex[j].x;
			edgeHalfVertex[out->sulci[i].edgeID[j]].y = out->sulci[i].vertex[j].y;
			edgeHalfVertex[out->sulci[i].edgeID[j]].z = out->sulci[i].vertex[j].z;
		}
	}

	geodesicThresh = 10.0;
	curvThresh = 0.0;

	newOut->numberOfSulci = out->numberOfSulci;
	newOut->sulci = (Sulci*)malloc(sizeof(Sulci)*(newOut->numberOfSulci+1));

	VtxOnFace start, end1, end2;

	for(i=0; i<out->numberOfSulci; i++)
	{
		if(!sulciStartVisited[i])
		{
			start.ID = -1;
			pathLen1 = 0;

			FindStart(surface,out->sulci[i].edgeID[0], out->sulci[i].vertex[0], &start);
			
			if(start.ID >= 0 && start.ID < surface->edgeNum)
			{
				end1.ID = -1;
			
				SearchConnectingHalfVertexUsingFastMarching(surface, segment, 
					start, marchingSpeed, edgeSulci, edgeHalfVertex, geodesicThresh, i, &end1);

				if( end1.ID >=0 && end1.ID < surface->edgeNum)
				{
					PathFinderHalfVertex(surface, marchingSpeed, start, end1, pathVtxList1, &pathLen1);

					//pathLen = 0;
					//AddSulciSegmentHalfVertex(surface, outTest, start, end, pathVtxList, pathLen);

					// check whether the path contains large curvature vertex  
					for(j=0; j<pathLen1; j++)
					{
						if(InterpolateCurvature(surface,pathVtxList1[j]) > curvThresh)
						{
							pathLen1 = 0;
							break;
						}
					}
					/*
					// the end1 is the end point or start point of another sulcus
					if(edgeSulciEnd[end1.ID] > 0)
					{
						if(edgeSulciEnd[end1.ID] == 1)
							sulciStartVisited[edgeSulci[end1.ID]] = true;
						else
							sulciEndVisited[edgeSulci[end1.ID]] = true;
					}
					*/
				}

			}
			sulciStartVisited[i] = true;
		}

		if(!sulciEndVisited[i])
		{
			start.ID = -1;
			pathLen2 = 0;

			FindStart(surface,out->sulci[i].edgeID[out->sulci[i].n-1], out->sulci[i].vertex[out->sulci[i].n-1], &start);
		
			if(start.ID >= 0 && start.ID < surface->edgeNum)
			{
				end2.ID = -1;
				SearchConnectingHalfVertexUsingFastMarching(surface, segment, 
					start, marchingSpeed, edgeSulci, edgeHalfVertex, geodesicThresh, i, &end2);

				if( end2.ID >=0 && end2.ID < surface->edgeNum)
				{
					PathFinderHalfVertex(surface, marchingSpeed, start, end2, pathVtxList2, &pathLen2);

					//pathLen = 0;
					//AddSulciSegmentHalfVertex(surface, outTest, start, end, pathVtxList, pathLen);

					// check whether the path contains large curvature vertex  
					for(j=0; j<pathLen2; j++)
					{
						if(InterpolateCurvature(surface,pathVtxList2[j]) > curvThresh)
						{
							pathLen2 = 0;
							break;
						}
					}
					/*					
					if(edgeSulciEnd[end2.ID] > 0)
					{
						if(edgeSulciEnd[end2.ID] == 1)
							sulciStartVisited[edgeSulci[end2.ID]] = true;
						else
							sulciEndVisited[edgeSulci[end2.ID]] = true;
					}
					*/
				}
			}
			sulciEndVisited[i] = true;
		}

		AddTwoSulciSegmentHalfVertex(surface, out, newOut, i, end1, pathVtxList1, pathLen1, end2, pathVtxList2, pathLen2);

		printf("i = %d\n",i);
	}
	
	free(edgeSulci);
	free(edgeSulciEnd);
	free(edgeHalfVertex);
	free(sulciStartVisited);
	free(sulciEndVisited);
	free(marchingSpeed);
	free(marchingSpeedNU);
	free(pathVtxList1);
	free(pathVtxList2);
	
	return;
}

void ConnectInterruptedSulciHalfVertexAndVertex(Surface *surface, Segment *segment, SulciTrackingOut *out, SulciTrackingOut *newOut)
{  
	int *edgeSulci;
	int *vertexSulci;

	int *edgeSulciOrder;
	int *vertexSulciOrder;

	unsigned char *edgeSulciEnd;

	Fvector3d *edgeHalfVertex;
	Fvector3d *vertexHalfVertex;

	unsigned char **connectSulci;
	unsigned char **connectSulciColor;
	int **edgeConnectSulci;
	
	double *marchingSpeed;
	double *marchingSpeedNU;
	
	VtxOnFace *pathVtxList1, *pathVtxList2;
	int pathLen1, pathLen2;
	int N;

	int i, j, k;
	float curvThresh;
	float geodesicThresh;
	
	// record sulci color in each edge
	edgeSulci = (int *)malloc(sizeof(int)*surface->edgeNum);

	// record sulci color in each vertex
	vertexSulci = (int *)malloc(sizeof(int)*surface->vertexNum);

	edgeSulciOrder = (int *)malloc(sizeof(int)*surface->edgeNum);
	vertexSulciOrder = (int *)malloc(sizeof(int)*surface->vertexNum);

	// record sulci position
	edgeHalfVertex = (Fvector3d *)malloc(sizeof(Fvector3d)*surface->edgeNum);
	vertexHalfVertex = (Fvector3d *)malloc(sizeof(Fvector3d)*surface->vertexNum);

	// recorde whether two sulci have been connected
	connectSulci = UCalloc2d(out->numberOfSulci, out->numberOfSulci);
	connectSulciColor = UCalloc2d(out->numberOfSulci, out->numberOfSulci);
	edgeConnectSulci = Ialloc2d(out->numberOfSulci, out->numberOfSulci);

	marchingSpeed = (double *)malloc(sizeof(double)*surface->vertexNum);
	marchingSpeedNU = (double *)malloc(sizeof(double)*surface->vertexNum);

	pathVtxList1 = (VtxOnFace *)malloc(sizeof(VtxOnFace)*1000);
	pathVtxList2 = (VtxOnFace *)malloc(sizeof(VtxOnFace)*1000);

	for(i=0; i<surface->vertexNum; i++)
	{
		vertexSulci[i] = -1;
		vertexSulciOrder[i] = -1;
		vertexHalfVertex[i].x = -1.0;
		vertexHalfVertex[i].y = -1.0;
		vertexHalfVertex[i].z = -1.0;
	}
	
	for(i=0; i<surface->edgeNum; i++)
	{
		edgeSulci[i] = -1;
		edgeSulciOrder[i] = -1;
		edgeHalfVertex[i].x = -1.0;
		edgeHalfVertex[i].y = -1.0;
		edgeHalfVertex[i].z = -1.0;
	}

	for(i=0; i<out->numberOfSulci; i++)
	for(j=0; j<out->numberOfSulci; j++)
	{
		connectSulci[i][j] = 0;
		connectSulciColor[i][j] = 0;
		edgeConnectSulci[i][j] = -1;
	}

	// uniform marching speed
	for(i=0; i<surface->vertexNum; i++)
	{
		marchingSpeed[i] = 1.0;

		if(surface->neighbors[i].size() != surface->adjacentfaces[i].size())
			marchingSpeed[i] = 0.00000001;
	}

	// nonuniform marching speed
	ComputeNonUniformlMarchingSpeed(surface, marchingSpeedNU);

	for(i=0; i<out->numberOfSulci; i++)
	{
		for(j=0; j<out->sulci[i].n; j++)
		{
			if(out->sulci[i].edgeID[j] >= 0)
			{
				edgeSulci[out->sulci[i].edgeID[j]] = i;
				edgeSulciOrder[out->sulci[i].edgeID[j]] = j;
				edgeHalfVertex[out->sulci[i].edgeID[j]].x = out->sulci[i].vertex[j].x;
				edgeHalfVertex[out->sulci[i].edgeID[j]].y = out->sulci[i].vertex[j].y;
				edgeHalfVertex[out->sulci[i].edgeID[j]].z = out->sulci[i].vertex[j].z;
			}
			else if(out->sulci[i].vertexID[j] >= 0)
			{
				vertexSulci[out->sulci[i].vertexID[j]] = i;
				vertexSulciOrder[out->sulci[i].vertexID[j]] = j;
				vertexHalfVertex[out->sulci[i].vertexID[j]].x = out->sulci[i].vertex[j].x;
				vertexHalfVertex[out->sulci[i].vertexID[j]].y = out->sulci[i].vertex[j].y;
				vertexHalfVertex[out->sulci[i].vertexID[j]].z = out->sulci[i].vertex[j].z;
			}
		}
	}

	curvThresh = 0;
	geodesicThresh = 10.0;

	newOut->numberOfSulci = out->numberOfSulci;
	newOut->sulci = (Sulci*)malloc(sizeof(Sulci)*(newOut->numberOfSulci+1));

	VtxOnFace start, end1, end2;

	for(i=0; i<out->numberOfSulci; i++)
	{
		pathLen1 = 0;
		end1.ID = -1;

		if(!out->sulci[i].isJunction[0] && (out->sulci[i].vertexID[0] >= 0 || out->sulci[i].edgeID[0] >= 0))
		{
			start.ID = -1;
			
			if(out->sulci[i].vertexID[0] >= 0)
			{
				start.lambda = 0.0;
				start.ID = out->sulci[i].vertexID[0];
			}
			else
				FindStart(surface,out->sulci[i].edgeID[0], out->sulci[i].vertex[0], &start);

			if(start.ID >= 0)
			{
				end1.ID = -1;
			
				SearchConnectingHalfVertexAndVertex(surface, segment, out, 
					start, marchingSpeed, edgeSulci, vertexSulci, edgeHalfVertex, geodesicThresh, out->sulci[i].color, &end1);

				if( end1.ID >=0)
				{
					bool endlessLoop = false;

					endlessLoop = PathFinderHalfVertex(surface, marchingSpeedNU, start, end1, pathVtxList1, &pathLen1);

					if(endlessLoop)
						pathLen1 = 0;

					if(end1.lambda == 0.0)
					{
						if((connectSulci[i][vertexSulci[end1.ID]] == 1 || connectSulciColor[out->sulci[i].color][out->sulci[vertexSulci[end1.ID]].color] == 1) )
							pathLen1 = 0;
					}
					else 
					{
						if(connectSulci[i][edgeSulci[end1.ID]] == 1|| connectSulciColor[out->sulci[i].color][out->sulci[edgeSulci[end1.ID]].color] == 1)
							pathLen1 = 0;
					}

					// check whether the path contains large curvature vertex  
					for(j=0; j<pathLen1; j++)
					{
						if(InterpolateCurvature(surface,pathVtxList1[j]) > curvThresh)
						{
							pathLen1 = 0;
							break;
						}

						if(pathVtxList1[j].lambda == 0.0)
						{
							if(( vertexSulci[pathVtxList1[j].ID] == i) || (j>0 &&  vertexSulci[pathVtxList1[j].ID] >=0))
							{
								pathLen1 = 0;
								break;
							}
						}
					}

					if(pathLen1 > 0)
					{
						if(end1.lambda == 0.0)
						{
							connectSulci[i][vertexSulci[end1.ID]] = 1;
							connectSulci[vertexSulci[end1.ID]][i] = 1;
		
							connectSulciColor[out->sulci[i].color][out->sulci[vertexSulci[end1.ID]].color] = 1;
							connectSulciColor[out->sulci[vertexSulci[end1.ID]].color][out->sulci[i].color] = 1;
							
							out->sulci[vertexSulci[end1.ID]].isJunction[vertexSulciOrder[end1.ID]] = true;
							out->sulci[vertexSulci[end1.ID]].junction = true;				
						}
						else
						{
							connectSulci[i][edgeSulci[end1.ID]] = 1;
							connectSulci[edgeSulci[end1.ID]][i] = 1;

							connectSulciColor[out->sulci[i].color][out->sulci[edgeSulci[end1.ID]].color] = 1;
							connectSulciColor[out->sulci[edgeSulci[end1.ID]].color][out->sulci[i].color] = 1;

							out->sulci[edgeSulci[end1.ID]].isJunction[edgeSulciOrder[end1.ID]] = true;
							out->sulci[edgeSulci[end1.ID]].junction = true;
						}
					}
				}
			}
		}

		pathLen2 = 0;
		end2.ID = -1;

		if(!out->sulci[i].isJunction[out->sulci[i].n-1] && (out->sulci[i].vertexID[out->sulci[i].n-1] >= 0 || out->sulci[i].edgeID[out->sulci[i].n-1] >= 0))
		{
			start.ID = -1;

			if(out->sulci[i].vertexID[out->sulci[i].n-1] >= 0)
			{
				start.lambda = 0.0;
				start.ID = out->sulci[i].vertexID[out->sulci[i].n-1];
			}
			else
                FindStart(surface,out->sulci[i].edgeID[out->sulci[i].n-1], out->sulci[i].vertex[out->sulci[i].n-1], &start);
		
			if(start.ID >= 0)
			{
				end2.ID = -1;
				SearchConnectingHalfVertexAndVertex(surface, segment, out, 
					start, marchingSpeed, edgeSulci, vertexSulci, edgeHalfVertex, geodesicThresh,out->sulci[i].color, &end2);

				if( end2.ID >=0 )
				{
					bool endlessLoop = false;

					endlessLoop = PathFinderHalfVertex(surface, marchingSpeedNU, start, end2, pathVtxList2, &pathLen2);

					if(endlessLoop)
						pathLen2 = 0;

					if( end2.lambda == 0.0 )
					{
						if( (connectSulci[i][vertexSulci[end2.ID]] == 1 || connectSulciColor[out->sulci[i].color][out->sulci[vertexSulci[end2.ID]].color]))
							pathLen2 = 0;
					}
					else 
					{
						if(connectSulci[i][edgeSulci[end2.ID]] == 1 || connectSulciColor[out->sulci[i].color][out->sulci[edgeSulci[end2.ID]].color] == 1)
							pathLen2 = 0;
					}
				
					for(j=0; j<pathLen2; j++)
					{
						if(InterpolateCurvature(surface,pathVtxList2[j]) > curvThresh)
						{
							pathLen2 = 0;
							break;
						}

						if(pathVtxList2[j].lambda == 0.0)
						{
							if( (vertexSulci[pathVtxList2[j].ID] == i) || (j>0 && vertexSulci[pathVtxList2[j].ID] >=0))
							{
								pathLen2 = 0;
								break;
							}
						}
					}

					if(pathLen2 > 0)
					{
						if(end2.lambda == 0.0)
						{
							connectSulci[i][vertexSulci[end2.ID]] = 1;
							connectSulci[vertexSulci[end2.ID]][i] = 1;
								
							connectSulciColor[out->sulci[i].color][out->sulci[vertexSulci[end2.ID]].color] = 1;
							connectSulciColor[out->sulci[vertexSulci[end2.ID]].color][out->sulci[i].color] = 1;

							out->sulci[vertexSulci[end2.ID]].isJunction[vertexSulciOrder[end2.ID]] = true;
							out->sulci[vertexSulci[end2.ID]].junction = true;
						}
						else
						{
							connectSulci[i][edgeSulci[end2.ID]] = 1;
							connectSulci[edgeSulci[end2.ID]][i] = 1;

							connectSulciColor[out->sulci[i].color][out->sulci[edgeSulci[end2.ID]].color] = 1;
							connectSulciColor[out->sulci[edgeSulci[end2.ID]].color][out->sulci[i].color] = 1;

							out->sulci[edgeSulci[end2.ID]].isJunction[edgeSulciOrder[end2.ID]] = true;
							out->sulci[edgeSulci[end2.ID]].junction = true;
						}
					}
				}
			}
		}

		AddTwoSulciSegmentHalfVertex(surface, out, newOut, i, end1, pathVtxList1, pathLen1, end2, pathVtxList2, pathLen2);

		printf("i = %d\n",i);
	}
	
	// reallocate color for temp sulci
	newOut->color2Sulci.resize(newOut->numberOfSulci);

	for(i=0; i<newOut->numberOfSulci; i++)
		newOut->color2Sulci[newOut->sulci[i].color].push_back(i);
	
	int count;

	for(i=0; i<newOut->numberOfSulci; i++)
	for(j=i+1; j<out->numberOfSulci; j++)
	{
		if(connectSulci[i][j] == 1&& newOut->sulci[i].color != newOut->sulci[j].color)
		{
			if(newOut->color2Sulci[newOut->sulci[i].color].size() >= newOut->color2Sulci[newOut->sulci[j].color].size())
			{
				count = newOut->color2Sulci[newOut->sulci[j].color].size();
				
				for(k=0; k<count; k++)
				{
					newOut->color2Sulci[newOut->sulci[i].color].push_back(newOut->color2Sulci[newOut->sulci[j].color][k]);
					newOut->color2Sulci[newOut->sulci[j].color].pop_back();
					newOut->sulci[newOut->color2Sulci[newOut->sulci[j].color][k]].color = newOut->sulci[i].color;
				}
			}
			else
			{
				count = newOut->color2Sulci[newOut->sulci[i].color].size();

				for(k=0; k<count; k++)
				{
					newOut->color2Sulci[newOut->sulci[j].color].push_back(newOut->color2Sulci[newOut->sulci[i].color][k]);
					newOut->color2Sulci[newOut->sulci[i].color].pop_back();
					newOut->sulci[newOut->color2Sulci[newOut->sulci[i].color][k]].color = newOut->sulci[j].color;					
				}
			}
		}
	}

	InterruptSulci(out, newOut);

	/*
	float length;
	int start, end;

	lenThresh = 8.0;

	newOut->numberOfSulci = 0;
	newOut->sulci = NULL;

	for(i=0; i<out->numberOfSulci; i++)
	{
		if(out->sulci[i].junction == 1)
		{
			start = 0;
			length = 0.0;
			
			for(j=start+1; j<out->sulci[i].n; j++)
			{
				if(out->sulci[i].isJunction[j] || j == out->sulci[i].n-1)
				{
					end = j;
					
					if(length > lenThresh)
						AddNewSulciFromOldSulci(out, newOut, i, start, end);
	
					start = end;
					length = 0;
				}
				length += DistanceBetweenTwoVertices(out->sulci[i].vertex[j-1], out->sulci[i].vertex[j]);
			}
		}
		else
			AddNewSulciFromOldSulci(out, newOut, i, 0, out->sulci[i].n-1);
	}
	*/

	// free memory
	free(edgeSulci);
	free(vertexSulci);
	
	free(edgeHalfVertex);
	free(vertexHalfVertex);
	
	free(marchingSpeed);
	free(marchingSpeedNU);
	
	free(pathVtxList1);
	free(pathVtxList2);
	
	Ifree2d(edgeConnectSulci,out->numberOfSulci);
	UCfree2d(connectSulci,out->numberOfSulci);
	UCfree2d(connectSulciColor,out->numberOfSulci);

	return;
}



void AddOneSulcus(Surface *surface, SulciTrackingOut *newOut, int sulcusLabel, int sulcusColor, VtxOnFace start, VtxOnFace end, VtxOnFace *pathVtxList, int pathLen)
{
	int i;
	Fvector3d V;

	newOut->sulci[sulcusLabel].n = pathLen + 1;

	newOut->sulci[sulcusLabel].color = sulcusColor;
	newOut->sulci[sulcusLabel].junction = false;
	newOut->sulci[sulcusLabel].curv1 = (float *)malloc(sizeof(float)*newOut->sulci[sulcusLabel].n);
	newOut->sulci[sulcusLabel].vertexID = (int *)malloc(sizeof(int)*newOut->sulci[sulcusLabel].n);
	newOut->sulci[sulcusLabel].edgeID = (int *)malloc(sizeof(int)*newOut->sulci[sulcusLabel].n);
	newOut->sulci[sulcusLabel].faceID = (int *)malloc(sizeof(int)*newOut->sulci[sulcusLabel].n);
	newOut->sulci[sulcusLabel].isJunction = (bool *)malloc(sizeof(bool)*newOut->sulci[sulcusLabel].n);
	newOut->sulci[sulcusLabel].vertex = (Fvector3d *)malloc(sizeof(Fvector3d)*newOut->sulci[sulcusLabel].n);

	if(pathLen > 0)
	{
		for(i=0; i<pathLen; i++)
		{
			if(pathVtxList[i].lambda == 0.0)
			{
				V.x = surface->vertex[pathVtxList[i].ID].x;
				V.y = surface->vertex[pathVtxList[i].ID].y;
				V.z = surface->vertex[pathVtxList[i].ID].z;
				
				newOut->sulci[sulcusLabel].vertexID[i] =  pathVtxList[i].ID;
				newOut->sulci[sulcusLabel].edgeID[i] = -1;
				newOut->sulci[sulcusLabel].faceID[i] = -1;
				//newOut->sulci[sulcusLabel].isJunction[i] = true;
				newOut->sulci[sulcusLabel].isJunction[i] = false;
			}
			else if(pathVtxList[i].lambda > 0.0)
			{
				V = InterpolateHalfVertex(surface->vertex[surface->edges[pathVtxList[i].ID].vertexID[0]], 
					surface->vertex[surface->edges[pathVtxList[i].ID].vertexID[1]], pathVtxList[i].lambda);
				
				newOut->sulci[sulcusLabel].vertexID[i] =  -1;
				newOut->sulci[sulcusLabel].edgeID[i] = pathVtxList[i].ID;
				newOut->sulci[sulcusLabel].faceID[i] =  -1;
				newOut->sulci[sulcusLabel].isJunction[i] = false;
			}
			else
			{
				V.x = (surface->vertex[surface->faces[pathVtxList[i].ID].v[0]].x 
					+ surface->vertex[surface->faces[pathVtxList[i].ID].v[1]].x + surface->vertex[surface->faces[pathVtxList[i].ID].v[2]].x)/3;
				V.y = (surface->vertex[surface->faces[pathVtxList[i].ID].v[0]].y 
					+ surface->vertex[surface->faces[pathVtxList[i].ID].v[1]].y + surface->vertex[surface->faces[pathVtxList[i].ID].v[2]].y)/3;
				V.z = (surface->vertex[surface->faces[pathVtxList[i].ID].v[0]].z 
					+ surface->vertex[surface->faces[pathVtxList[i].ID].v[1]].z + surface->vertex[surface->faces[pathVtxList[i].ID].v[2]].z)/3;

				newOut->sulci[sulcusLabel].vertexID[i] =  -1;
				newOut->sulci[sulcusLabel].edgeID[i] = -1;
				newOut->sulci[sulcusLabel].faceID[i] = pathVtxList[i].ID;
				newOut->sulci[sulcusLabel].isJunction[i] = false;
			}
			newOut->sulci[sulcusLabel].curv1[i] =  0.0;
			newOut->sulci[sulcusLabel].vertex[i].x = V.x;
			newOut->sulci[sulcusLabel].vertex[i].y = V.y;
			newOut->sulci[sulcusLabel].vertex[i].z = V.z;
		}

		if(end.lambda == 0.0)
		{
			V.x = surface->vertex[end.ID].x;
			V.y = surface->vertex[end.ID].y;
			V.z = surface->vertex[end.ID].z;
			
			newOut->sulci[sulcusLabel].vertexID[pathLen] =  end.ID;
			newOut->sulci[sulcusLabel].edgeID[pathLen] = -1;
			newOut->sulci[sulcusLabel].faceID[pathLen] = -1;
			//newOut->sulci[sulcusLabel].isJunction[pathLen] = true;
			newOut->sulci[sulcusLabel].isJunction[pathLen] = false;
		}			
		else if(end.lambda > 0.0)
		{
			V = InterpolateHalfVertex(surface->vertex[surface->edges[end.ID].vertexID[0]], 
			surface->vertex[surface->edges[end.ID].vertexID[1]], end.lambda);

			newOut->sulci[sulcusLabel].vertexID[pathLen] =  -1;
			newOut->sulci[sulcusLabel].edgeID[pathLen] = end.ID;
			newOut->sulci[sulcusLabel].faceID[pathLen] = -1;
			newOut->sulci[sulcusLabel].isJunction[pathLen] = false;
		}
		else
		{
			V.x = (surface->vertex[surface->faces[end.ID].v[0]].x 
				+ surface->vertex[surface->faces[end.ID].v[1]].x + surface->vertex[surface->faces[end.ID].v[2]].x)/3;
			V.y = (surface->vertex[surface->faces[end.ID].v[0]].y 
				+ surface->vertex[surface->faces[end.ID].v[1]].y + surface->vertex[surface->faces[end.ID].v[2]].y)/3;
			V.z = (surface->vertex[surface->faces[end.ID].v[0]].z 
				+ surface->vertex[surface->faces[end.ID].v[1]].z + surface->vertex[surface->faces[end.ID].v[2]].z)/3;
			
			newOut->sulci[sulcusLabel].vertexID[pathLen] =  -1;
			newOut->sulci[sulcusLabel].edgeID[pathLen] = -1;
			newOut->sulci[sulcusLabel].faceID[pathLen] = end.ID;
			newOut->sulci[sulcusLabel].isJunction[pathLen] = false;
		}
		newOut->sulci[sulcusLabel].curv1[pathLen] = 0.0;
		newOut->sulci[sulcusLabel].vertex[pathLen].x = V.x;
		newOut->sulci[sulcusLabel].vertex[pathLen].y = V.y;
		newOut->sulci[sulcusLabel].vertex[pathLen].z = V.z;
	}

	return;
}


void KeepOneSulcus(Surface *surface, SulciTrackingOut *newOut, SulciTrackingOut *out, int sulcusLabel, int sulcusColor)
{
	int i;

	newOut->sulci[sulcusLabel].n = out->sulci[sulcusLabel].n;
	newOut->sulci[sulcusLabel].junction = false;
	newOut->sulci[sulcusLabel].color = out->sulci[sulcusLabel].color;
	newOut->sulci[sulcusLabel].curv1 = (float *)malloc(sizeof(float)*newOut->sulci[sulcusLabel].n);
	newOut->sulci[sulcusLabel].vertexID = (int *)malloc(sizeof(int)*newOut->sulci[sulcusLabel].n);
	newOut->sulci[sulcusLabel].edgeID = (int *)malloc(sizeof(int)*newOut->sulci[sulcusLabel].n);
	newOut->sulci[sulcusLabel].faceID = (int *)malloc(sizeof(int)*newOut->sulci[sulcusLabel].n);
	newOut->sulci[sulcusLabel].isJunction = (bool *)malloc(sizeof(bool)*newOut->sulci[sulcusLabel].n);
	newOut->sulci[sulcusLabel].vertex = (Fvector3d *)malloc(sizeof(Fvector3d)*newOut->sulci[sulcusLabel].n);

	for(i=0; i<newOut->sulci[sulcusLabel].n; i++)
	{
		newOut->sulci[sulcusLabel].curv1[i] = out->sulci[sulcusLabel].curv1[i];
		newOut->sulci[sulcusLabel].vertex[i].x = out->sulci[sulcusLabel].vertex[i].x ;
		newOut->sulci[sulcusLabel].vertex[i].y = out->sulci[sulcusLabel].vertex[i].y ;
		newOut->sulci[sulcusLabel].vertex[i].z = out->sulci[sulcusLabel].vertex[i].z ;
		newOut->sulci[sulcusLabel].vertexID[i] = out->sulci[sulcusLabel].vertexID[i];
		newOut->sulci[sulcusLabel].edgeID[i] = out->sulci[sulcusLabel].edgeID[i];
		newOut->sulci[sulcusLabel].faceID[i] = newOut->sulci[sulcusLabel].faceID[i];
		newOut->sulci[sulcusLabel].isJunction[i] = false;
	}

	return;
}



void ReconnectStartEndUsingFastMarching(Surface *surface, Segment *segment, SulciTrackingOut *out, SulciTrackingOut *newOut)
{  
	int N;
	int i, j;
	int pathLen;
	float geodesicThresh;
	
	double *marchingSpeedNU;
	VtxOnFace *pathVtxList;
	VtxOnFace start, end;

	marchingSpeedNU = (double *)malloc(sizeof(double)*surface->vertexNum);
	pathVtxList = (VtxOnFace *)malloc(sizeof(VtxOnFace)*2000);

	newOut->numberOfSulci = out->numberOfSulci;
	newOut->sulci = (Sulci*)malloc(sizeof(Sulci)*(newOut->numberOfSulci+1));

	ComputeNonUniformlMarchingSpeedNew(surface, segment, marchingSpeedNU);

	for(i=0; i<out->numberOfSulci; i++)
	{
		if(out->sulci[i].vertexID[0] >= 0)
		{
			start.lambda = 0.0;
			start.ID = out->sulci[i].vertexID[0];
		}
		else if(out->sulci[i].edgeID[0] >= 0)
			FindStart(surface,out->sulci[i].edgeID[0], out->sulci[i].vertex[0], &start);
		else
			FindStart(surface,out->sulci[i].edgeID[1], out->sulci[i].vertex[1], &start);

		if(out->sulci[i].vertexID[out->sulci[i].n-1] >= 0)
		{
			end.lambda = 0.0;
			end.ID = out->sulci[i].vertexID[out->sulci[i].n-1] ;
		}
		else if(out->sulci[i].edgeID[out->sulci[i].n-1] >= 0)
			FindStart(surface,out->sulci[i].edgeID[out->sulci[i].n-1], out->sulci[i].vertex[out->sulci[i].n-1], &end);
		else
			FindStart(surface,out->sulci[i].edgeID[out->sulci[i].n-2], out->sulci[i].vertex[out->sulci[i].n-2], &end);

		bool endlessLoop = false;
		pathLen = 0;
		endlessLoop = PathFinderHalfVertex(surface, marchingSpeedNU, end, start, pathVtxList, &pathLen);
		
		if(out->sulci[i].vertexID[0] < 0 && out->sulci[i].edgeID[0] < 0)
		{
			for(j=pathLen-1; j>=0; j--)
			{
				pathVtxList[j+1].ID = pathVtxList[j].ID;
				pathVtxList[j+1].lambda = pathVtxList[j].lambda;
			}

			pathVtxList[0].ID = out->sulci[i].faceID[0];
			pathVtxList[0].lambda = -1.0;

			pathLen += 1;
		}
		if(out->sulci[i].vertexID[out->sulci[i].n-1] < 0 && out->sulci[i].edgeID[out->sulci[i].n-1] < 0)
		{
			pathVtxList[pathLen].ID = end.ID;
			pathVtxList[pathLen].lambda = end.lambda;

			end.lambda = -1.0;
			end.ID = out->sulci[i].faceID[out->sulci[i].n-1];
			
			pathLen += 1;
		}

		
		if(!endlessLoop)
			AddOneSulcus(surface, newOut, i, out->sulci[i].color, start, end, pathVtxList, pathLen);
		else
			KeepOneSulcus(surface, newOut, out, i, out->sulci[i].color);

		printf("i = %d\n",i);
	}

	free(marchingSpeedNU);
	free(pathVtxList);
	
	return;
}


void ChangePincinpalDirection(Surface *surface)
{
	int i;
	
	for(i=0; i<surface->vertexNum; i++)
	{
		if(surface->dcurv[0][i] > 0)
		{
			surface->pdir1[i].x = -1.0*surface->pdir1[i].x;
			surface->pdir1[i].y = -1.0*surface->pdir1[i].y;
			surface->pdir1[i].z = -1.0*surface->pdir1[i].z;
			surface->dcurv[0][i] = -1.0*surface->dcurv[0][i];
		}
	}

	return;
}

void SwapInt(int *x, int *y)
{
	int temp;
	
	temp = *y; *y = *x; *x = temp;

	return;
}

int *CreateLookupTable(int numberOfSulci)
{
	int i;
	int *lookupTable;

	lookupTable = (int *)malloc(sizeof(int)*(numberOfSulci));

	for(i=0; i<numberOfSulci; i++)
		lookupTable[i] = i;

	 srand ( time(NULL) );

	for(i=0; i<numberOfSulci; i++)
		SwapInt(&lookupTable[i],&lookupTable[(int)(rand()*(numberOfSulci)/(RAND_MAX+1.0))]);

	return  lookupTable;
}

void WriteSulciVTK(char *fname, SulciTrackingOut *out, int *lookupTable)
{
	FILE *fp;

	int i, t, pointNum;
	int position;

	Fvector3d vector;

	pointNum = 0;
	for(t=0; t<out->numberOfSulci; t++)
		pointNum += out->sulci[t].n;

	fp = fopen(fname,"wt");

	fprintf(fp,"# vtk DataFile Version 3.0\n");
	fprintf(fp,"Sulci\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp,"POINTS %d float\n",pointNum);
	
	for(t=0; t<out->numberOfSulci; t++)
	{
		for(i=0; i<out->sulci[t].n; i++)
		{
			fprintf(fp,"%f %f %f\n",out->sulci[t].vertex[i].x, out->sulci[t].vertex[i].y, out->sulci[t].vertex[i].z);
		}
	}

	position = 0;
	fprintf(fp,"CELLS %d %d\n", out->numberOfSulci, pointNum+out->numberOfSulci);
	for(t=0; t<out->numberOfSulci; t++)
	{
		fprintf(fp,"%d ", out->sulci[t].n);

		for(i=0; i<out->sulci[t].n; i++)
		{
			fprintf(fp,"%d ",position+i);
		}
		fprintf(fp,"\n");

		position += out->sulci[t].n;
	}

	fprintf(fp,"CELL_TYPES %d\n",out->numberOfSulci);
	for(t=0; t<out->numberOfSulci; t++)
	{
		fprintf(fp,"4\n");
	}

	// curv1 as lookup table
	fprintf(fp,"POINT_DATA %d\n",pointNum);
	fprintf(fp,"SCALARS scalars float\n");
	fprintf(fp,"LOOKUP_TABLE default\n");
	
	for(t=0; t<out->numberOfSulci; t++)
	{
		for(i=0; i<out->sulci[t].n; i++)
		{
			fprintf(fp,"%f\n",out->sulci[t].curv1[i]);
		}
	}

	// sulci length as lookup table
	float *length;
	
	length = (float *)malloc(sizeof(float)*out->numberOfSulci);

	for(t=0; t<out->numberOfSulci; t++)
		length[t] = 0.0;

	for(t=0; t<out->numberOfSulci; t++)
	{
		for(i=0; i<out->sulci[t].n-1; i++)
			length[t] += DistanceBetweenTwoVertices(out->sulci[t].vertex[i], out->sulci[t].vertex[i+1]);
	}

	fprintf(fp,"SCALARS length float\n");
	fprintf(fp,"LOOKUP_TABLE lengthTable\n");
	for(t=0; t<out->numberOfSulci; t++)
	{
		for(i=0; i<out->sulci[t].n; i++)
		{
			fprintf(fp,"%f\n",length[t]);
		}
	}
	free(length);

	// random color code sulci

	fprintf(fp,"SCALARS color1 int\n");
	fprintf(fp,"LOOKUP_TABLE colorTable\n");
	for(t=0; t<out->numberOfSulci; t++)
	{
		for(i=0; i<out->sulci[t].n; i++)
		{
			fprintf(fp,"%d\n",lookupTable[t]);
		}
	}

	fprintf(fp,"SCALARS color2 int\n");
	fprintf(fp,"LOOKUP_TABLE colorTable\n");
	for(t=0; t<out->numberOfSulci; t++)
	{
		for(i=0; i<out->sulci[t].n; i++)
		{
			fprintf(fp,"%d\n",lookupTable[out->sulci[t].color]);
		}
	}

	fprintf(fp,"SCALARS edgeID int\n");
	fprintf(fp,"LOOKUP_TABLE edgeIDTable\n");
	for(t=0; t<out->numberOfSulci; t++)
	{
		for(i=0; i<out->sulci[t].n; i++)
		{
			fprintf(fp,"%d\n",out->sulci[t].edgeID[i]);
		}
	}

	/*
	fprintf(fp,"SCALARS isJunction int\n");
	fprintf(fp,"LOOKUP_TABLE isJunctionTable\n");
	for(t=0; t<out->numberOfSulci; t++)
	{
		for(i=0; i<out->sulci[t].n; i++)
		{
			fprintf(fp,"%d\n",out->sulci[t].isJunction[i]);
		}
	}
	*/
	/*
	fprintf(fp,"VECTORS pdir2 float\n");
	for(t=0; t<out->numberOfSulci; t++)
	{
		for(i=0; i<out->sulci[t].n; i++)
		{
			vector.x = out->sulci[t].pdir2[i].x;
			vector.y = out->sulci[t].pdir2[i].y;
			vector.z = out->sulci[t].pdir2[i].z;
			
			fprintf(fp,"%f %f %f\n",vector.x, vector.y, vector.z);
		}		
	}
*/
	fclose(fp);

	return;
}


void RealseSuface(Surface *surface)
{
	free(surface->vertex);
	free(surface->faces);
	free(surface->curv1);
	free(surface->curv2);
	//free(surface->depth);
	free(surface->pdir1);
	free(surface->pdir2);
	free(surface->normal);
	free(surface->isSulci);
	free(surface->across_edge);

	return;
}


int main(int argc, char *argv[])
{
	char *fname;
	char fnameSulci[200];
	char fnameSulciPrune[200];
	char fnameSulciPruneConnect[200];
	char fnameSulciPruneReconnect[200];

	SulciTrackingOut out;

	Params params;

	Segment segment;

	int *lookupTable;

	if(argc < 2)
	{
		printf("Usage: SulciTracking input.vtk\n");

		exit(1);
	}

	fname = argv[1];

	out.numberOfSulci = 0;

	ReadSurfaceAttribute(fname, &(out.surface));

	ChangePincinpalDirection(&(out.surface));

	GetSurfaceEdgesAttribute(&(out.surface));

	ComputeSulciSegment(&(out.surface), &params, &segment, false);

	ConnectSulciFromSegment(&out, &segment);
	
	RemoveSingleCenterPoint(&out);

	RemoveCandiateStartingSegment(&out, &segment);

//	sprintf(fnameSulci,"%s.StartSulci.vtk",fname);
//	WriteSegmentVTK(fnameSulci, &segment);

	lookupTable = CreateLookupTable(out.numberOfSulci);

	InitializeColor2SulciTable(&out);

//	sprintf(fnameSulciPrune,"%s.SulciPrune.vtk",fname);
//	WriteSulciVTK(fnameSulciPrune, &out, lookupTable);

	free(lookupTable);
	
	SulciTrackingOut outCombine;

	CombineAdjacentSulci(&out, &segment, &outCombine);

	lookupTable = CreateLookupTable(outCombine.numberOfSulci+100);

	sprintf(fnameSulciPruneConnect,"%s.SulciPruneCombine.vtk",fname);
	WriteSulciVTK(fnameSulciPruneConnect, &outCombine, lookupTable);

	SulciTrackingOut outReconnect;
	ReconnectStartEndUsingFastMarching(&(out.surface), &segment, &outCombine, &outReconnect);
	
	sprintf(fnameSulciPruneReconnect,"%s.SulciPruneCombineReconnect.vtk",fname);
	WriteSulciVTK(fnameSulciPruneReconnect, &outReconnect, lookupTable);

	SulciTrackingOut outTest;
	ConnectInterruptedSulciHalfVertexAndVertex(&(out.surface), &segment,  &outReconnect, &outTest);

	sprintf(fnameSulciPruneConnect,"%s.SulciPruneCombineReconnectConnect.vtk",fname);
	WriteSulciVTK(fnameSulciPruneConnect, &outTest, lookupTable);

	free(lookupTable);

	return 0;
}
