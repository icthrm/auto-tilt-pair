#include "ifgt.h"

// #define PT2DISTANCE(a,b)  (((a).x-(b).x)*((a).x-(b).x)+((a).y-(b).y)*((a).y-(b).y))

inline static float PT2DISTANCE(const util::point2d& a, const util::point2d& b)
{
	float deltax = ((a).x-(b).x);
	float deltay = ((a).y-(b).y);
	return deltax*deltax+deltay*deltay;
}

inline static int nchoosek(int n, int k)
{
    int k_orig = k;
    int n_k = n-k;
    if(k < n_k){
        k = n_k;
        n_k = n-k;
    }
    
    double nchsk = 1;
    for(int i = 1; i <= n_k; i++){
        nchsk *= ++k;
        nchsk /= i;
    }
    
    if(nchsk > __INT_MAX__){
        std::cout<<"n choose k for " <<n<< " and " << k_orig
           << " caused an overflow. Dimensionality of the data might be "
              "too high.";
		exit(0);
    }
    return int(nchsk);
}

inline static void choose_ifgt_params(double bandwidth, double epsilon, int max_num_clusters, int truncation_number_ul, int* nclusters, double* cutoff_radius)
{
	int dim = 2;
    double h2 = bandwidth*bandwidth;
    *cutoff_radius = bandwidth*std::sqrt(std::log(1.0/epsilon));
    double complexity_min = __DBL_MAX__;
    *nclusters = 0;

    for(int i = 0; i < max_num_clusters; i++){
        double rx = std::pow(double(i + 1), -1.0/double(dim));
        double rx2 = rx*rx;
        double n = std::min(double(i + 1), std::pow(*cutoff_radius/rx, dim));
        double error = __DBL_MAX__;
        double temp = 1.0;
        int p = 0;
        while ((error > epsilon) && (p <= truncation_number_ul)){
            ++p;
            double b = std::min((rx+std::sqrt(rx2+2.0*p*h2))*0.5, rx+*cutoff_radius);		//*0.5 = /2
            double c = rx-b;
            temp *= 2*rx*b/h2/p;
            error = temp*std::exp(-c*c/h2);
        }
        double complexity = i+1+std::log(double(i + 1))+(n+1)*nchoosek(p - 1 + dim, dim);
        if(complexity < complexity_min){
            complexity_min = complexity;
            *nclusters = i + 1;
        }
    }
}

// Chooses the appropriate truncation number for IFGT, given a max clustering radius.
inline static int choose_ifgt_truncation_number(double bandwidth, double epsilon, double rx, int truncation_number_ul)
{
	int dim = 2;
    double h2 = bandwidth*bandwidth;
    double rx2 = rx*rx;
    double r = std::min(std::sqrt(dim), bandwidth*std::sqrt(std::log(1.0/epsilon)));
    double error = __DBL_MAX__;
    int p = 0;
    double temp = 1.0;
	
    while((error > epsilon) && (p <= truncation_number_ul)){
        ++p;
        double b = std::min((rx+std::sqrt(rx2+2*double(p)*h2))*0.5, rx+r);
        double c = rx-b;
        temp *= 2*rx*b/h2/double(p);
        error = temp*std::exp(-(c*c)/h2);
    }
    
    return p;
}

void ifgt::IFGaussianTransform::CalculateConstantSeries(std::vector<double>& monomials)
{
	int dim = 2;
    std::vector<size_t> heads(dim+1, 0);
    heads[dim] = __INT_MAX__;
    std::vector<size_t> cinds(p_max_total, 0);
    monomials = std::vector<double>(p_max_total, 1);

    for(size_t k = 1, t = 1, tail = 1; k < truncation_number; k++, tail = t){
        for(size_t i = 0; i < dim; i++){
            size_t head = heads[i];
            heads[i] = t;
            for(size_t j = head; j < tail; ++j, ++t){
                cinds[t] = (j < heads[i + 1]) ? cinds[j] + 1 : 1;
                monomials[t] = 2.0*monomials[j];
                monomials[t] /= double(cinds[t]);
            }
        }
    }
}

void ifgt::IFGaussianTransform::CalculateMonomials(const double* d, std::vector<double>& monomials) const
{
	int dim = 2;
    std::vector<int> heads(dim, 0);
    new (&monomials) std::vector<double>(p_max_total, 1);
    for(int k = 1, t = 1, tail = 1; k < truncation_number; k++, tail = t){
        for(int i = 0; i < dim; ++i){
            int head = heads[i];
            heads[i] = t;
            for(int j = head; j < tail; ++j, ++t){
                monomials[t] = d[i]*monomials[j];
            }
        }
    }
}

ifgt::IFGaussianTransform::IFGaussianTransform(const std::vector<util::point2d>& source, double bw, double eps): source(source),
						bandwidth(bw), epsilon(eps), nclusters(0), truncation_number(0), p_max_total(0), constant_series()
{
	int max_num_clusters = round(0.2 * 100 / bandwidth);

	choose_ifgt_params(bandwidth, eps, max_num_clusters, TRUNCATION_NUMBER_UL, &nclusters, &cutoff_radius);
	
	std::cout<<nclusters<<std::endl;
    if (nclusters == 0) {
		std::cout<<"Don't need IFGT technique."<<std::endl;
    }
    
	new (&clustering)Clustering(source, nclusters, epsilon);			
    truncation_number = choose_ifgt_truncation_number(bandwidth, epsilon, clustering.max_radius, TRUNCATION_NUMBER_UL);
    p_max_total = nchoosek(truncation_number-1+2, 2);
    CalculateConstantSeries(constant_series);
	
    ry_square.resize(nclusters);
    for(int j = nclusters; j--;){		//for(int j = 0; j < nclusters; j++){
        double ry = cutoff_radius + clustering.radii[j];
        ry_square[j] = ry*ry;
    }
}

void ifgt::IFGaussianTransform::CalculateTransformWith(const std::vector< util::point2d >& target, std::vector<double>& result)
{
	std::vector<double> weights = std::vector<double>(source.size(), 1);
	return CalculateTransformWith(target, weights, result);
}

void ifgt::IFGaussianTransform::CalculateTransformWith(const std::vector<util::point2d>& target, const std::vector<double>& weights, std::vector<double>& result)
{
    double h2 = bandwidth*bandwidth;
	double h2_1 = 1/h2;
	
	double C[nclusters][p_max_total];
	memset(C, 0, sizeof(double)*nclusters*p_max_total);
	
    for(size_t i = source.size(); i--;){		//for(size_t i = 0; i < fixed.size(); i++){
        double distance = 0.0;
        double dx[2];
		
		double deltax = source[i].x-clustering.centers[clustering.indices[i]].x;
		double deltay = source[i].y-clustering.centers[clustering.indices[i]].y;
		distance = deltax*deltax+deltay*deltay;
		dx[0] = deltax/bandwidth;
		dx[1] = deltay/bandwidth;
		
        std::vector<double> monomials;
		CalculateMonomials(dx, monomials);
		
        double f = weights[i]*std::exp(-distance*h2_1);			//double f = weights[i]*std::exp(-distance/h2);
        for(int alpha = p_max_total; alpha--;){		//for(int alpha = 0; alpha < p_max_total; alpha++){
            C[clustering.indices[i]][alpha] += f*monomials[alpha];
        }
    }

    for(int j = nclusters; j--;){			//for(int j = 0; j < nclusters; j++){
        for(int alpha = p_max_total; alpha--;){//for(int alpha = 0; alpha < p_max_total; alpha++){
            C[j][alpha] *= constant_series[alpha];
        }
    }
	
    new(&result)std::vector<double>(target.size(), 0);

    for(size_t i = target.size(); i--;){		//for(size_t i = 0; i < target.size(); i++){
        for(int j = nclusters; j--;){		//for(int j = 0; j < nclusters; j++){
            double distance = 0.0;
            double dy[2];
			
			double deltax = target[i].x-clustering.centers[j].x;
			double deltay = target[i].y-clustering.centers[j].y;
			distance = deltax*deltax+deltay*deltay;
			if(distance <= ry_square[j]){
				dy[0] = deltax/bandwidth;
				dy[1] = deltay/bandwidth;
				std::vector<double> monomials;
				CalculateMonomials(dy, monomials);
				
				double g = std::exp(-distance*h2_1);		//std::exp(-distance/h2);
				double acc = 0;
				for(int k = monomials.size(); k--;){		//for(int k = 0; k < monomials.size(); k++){
					acc += C[j][k]*monomials[k];
				}
				result[i] += acc*g;
			}
        }
    }
}

static std::vector<util::point2d> generate_centers(const std::vector<util::point2d>& points, int nclusters)
{
	size_t size = points.size();
	std::vector<util::point2d> centers;
	int seed[nclusters];
	seed[0] = rand()%size;
	int idx = 1;
	
	while(true){
pickcluster:
		int tmp = rand()%size;
		for(int i = 0; i < idx; i++){
			if(tmp == seed[i]){
				goto pickcluster;
			}
		}
		seed[idx] = tmp;
		if(++idx >= nclusters){
			break;
		}
	}
	
	for(int i = 0; i < nclusters; i++){
		centers.push_back(points[seed[i]]);
	}
	
    return centers;
}

ifgt::Clustering::Clustering(const std::vector<util::point2d>& points, int nclusters, double epsilon)
{
    std::vector<util::point2d> icenters = generate_centers(points, nclusters);
    new (this)Clustering(points, epsilon, icenters);
}

ifgt::Clustering::Clustering(const std::vector<util::point2d>& points, double epsilon, const std::vector<util::point2d>& icenters)
{
	int nclusters = icenters.size();
	indices.resize(points.size());
	centers = icenters;
	double error = 0.0;
    double old_error = 0.0;
	
	do{
		npoints = std::vector<int>(nclusters, 0);
		std::vector<util::point2d> local_centers(centers.size(), util::point2d(0,0));
		old_error = error;
		error = 0;
		
		for(size_t i = points.size(); i--;){		//for(size_t i = 0; i < points.size(); i++){
			double min_distance = __DBL_MAX__;
			for(size_t j = centers.size(); j--;){		//for(size_t j = 0; j < centers.size(); j++){
				double distance = PT2DISTANCE(points[i], centers[j]);
				if(distance < min_distance){
					indices[i] = j;
					min_distance = distance;
				}
			}

			local_centers[indices[i]].x += points[i].x;
			local_centers[indices[i]].y += points[i].y;
			++npoints[indices[i]];
			error += min_distance;
		}
		
		for(size_t j = nclusters; j--;){		//for(size_t j = 0; j < nclusters; j++){
			local_centers[j].x /= npoints[j];
			local_centers[j].y /= npoints[j];
		}
		centers = local_centers;
	}while(std::abs(error-old_error) > epsilon);
	
	max_radius = __DBL_MIN__;
	radii = std::vector<double>(nclusters, __DBL_MIN__);
	
	for(size_t i = 0; i < points.size(); i++){
		double distance = sqrt(PT2DISTANCE(points[i],centers[indices[i]]));
		if(distance > radii[indices[i]]){
            radii[indices[i]] = distance;
        }
        if(distance > max_radius){
            max_radius = distance;
        }
	}
}
