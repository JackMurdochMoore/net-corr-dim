% Load an empirical network from data.
% 
% Associated with 
% "Correlation dimension in empirical networks" 
% by 
% Jack Murdoch Moore, Haiying Wang, Michael Small, Gang Yan, Huijie Yang, 
% and Changgui Gu. 
%
function [A, nameStr, nameStr2, networkFile] = load_network(networkFlag)

nameStr2 = '';

networkFolder = 'networks';

switch networkFlag
    
    case 1
        
        % % Attractive result:
        %
        % TV show, Network Data Repository
        % https://networkrepository.com/fb-pages-tvshow.php
        % N = 3,892
        % 
        % Rozemberczki (2019) GEMSEC: Graph Embedding with Self Clustering
        % From
        % https://networkrepository.com/fb-pages-tvshow.php
        % "Mutually liked facebook pages. Nodes represent the pages and
        % edges are mutual likes among them." 
        %
        networkFile = 'fb-pages-tvshow.edges';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d,%d', 'HeaderLines', 0); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency; N = size(A, 1); A(eye(N) == 1) = 0; G = graph(A);
        nameStr = 'FB pages, TV show';
        nameStr2 = 'Facebook mutual likes: TV show';
        
    case 2
        
        % % Attractive result:
        %
        % Power grid from Network Science, by Barabasi et al.
        % http://networksciencebook.com/translations/en/resources/data.html
        % N = 4,941
        %
        % Watts and Strogatz (1998) Collective dynamics of 'small-world'
        % networks:
        % "the electrical power grid of the western United States,"
        % 
        networkFile = 'powergrid.edgelist.txt';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d'); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency; N = size(A, 1); A(eye(N) == 1) = 0; G = graph(A);
        nameStr = 'Network science, power grid';
        nameStr2 = 'Electrical power grid';
        
    case 3
        
        % % Attractive result:
        %
        % Politician, Network Data Repository
        % https://networkrepository.com/fb-pages-politician.php
        % N = 5,908
        % 
        % Rozemberczki (2019) GEMSEC: Graph Embedding with Self Clustering
        % From
        % https://networkrepository.com/fb-pages-politician.php
        % "Mutually liked facebook pages. Nodes represent the pages and
        % edges are mutual likes among them." 
        %
        networkFile = 'fb-pages-politician.edges';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d,%d', 'HeaderLines', 0); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency; N = size(A, 1); A(eye(N) == 1) = 0; G = graph(A);
        nameStr = 'FB pages, politician';
        nameStr2 = 'Facebook mutual likes: Politician';
        
    case 12
        
        % % Attractive result
        % 
        % PhDs in computer science
        % https://networkrepository.com/ca-CSphd.php
        % N = 1025
        %
        % Batagelj et al. (2005) Exploratory social network analysis with
        % Pajek, p. 252:
        % "The file PhD.net contains the ties between Ph.D. students and
        % their advisors in theoretical computer science; each arc points
        % from an advisor to a student."  
        %
        networkFile = 'ca-CSphd.mtx';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d', 'HeaderLines', 3); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency;
        nameStr = 'CS PhD';
        nameStr2 = 'PhD supervision: Computer science';
    
    case 13
        
        % % Attractive result
        % 
        % Erdos collaboration network
        % https://networkrepository.com/ca-Erdos992.php
        % N = 4991
        %
        % Batagelj and Mrvar (2000) Some analyses of Erdos collaboration
        % graph:
        % "By removing Paul Erdos himself and connections to him from the
        % graph EE we get the truncated Erdos collaboration graph EEX. The
        % last, 1999 edition of this graph contains 6100 vertices and 9939
        % edges."
        %
        networkFile = 'ca-Erdos992.mtx';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d', 'HeaderLines', 3); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency;
        nameStr = 'Erdos';
        nameStr2 = 'Collaboration: Erdos';
        
    case 15
        
        % % Attractive result
        % 
        % Collaboration network of arXiv General Relativity
        % https://networkrepository.com/ca-GrQc.php
        % N = 4,158
        % 
        % J. Leskovec, J. Kleinberg and C. Faloutsos (2003) Graph
        % Evolution: Densification and Shrinking Diameters.
        % From
        % https://snap.stanford.edu/data/ca-GrGc.html
        % "Arxiv GR-QC (General Relativity and Quantum Cosmology)
        % collaboration network is from the e-print arXiv and covers
        % scientific collaborations between authors papers submitted to
        % General Relativity and Quantum Cosmology category. If an author i
        % co-authored a paper with author j, the graph contains a
        % undirected edge from i to j. If the paper is co-authored by k
        % authors this generates a completely connected (sub)graph on k
        % nodes. The data covers papers in the period from January 1993 to
        % April 2003 (124 months). It begins within a few months of the
        % inception of the arXiv, and thus represents essentially the
        % complete history of its GR-QC section." 
        %
        networkFile = 'ca-GrQc.mtx';
        fileID = fopen([networkFolder, '/', networkFile]);
        C = textscan(fileID, '%d %d', 'HeaderLines', 2);
        fclose(fileID);
        from_list = C{:, 1}; to_list = C{:, 2};
        minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency;
        nameStr = 'arXiv GR';
        nameStr2 = 'Collaboration: General relativity';
        
    case 16
        
        % % Attractive result:
        %
        % Udry schools data #27, largest connected component
        % https://github.com/JackMurdochMoore/small-world
        % N = 1,152
        % 
        % Bearman et al. (2004) Chains of affection: the structure of
        % adolescent romantic and sexual networks
        % 
        networkFile = 'udry_adol_health_large_conn_comp.mat';
        load([networkFolder, '/', networkFile], 'GCell')
        iiG = 27; G = GCell{iiG};
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency;
        nameStr = 'School 27';
        nameStr2 = 'High school friendship';
        
    case 17
        
        % % Attractive result:
        % 
        % Udry schools data #67, largest connected component
        % https://github.com/JackMurdochMoore/small-world
        % https://web.archive.org/web/20210115130726/http://moreno.ss.uci.edu/data.html#adhealth
        % N = 439
        % 
        % It is hard to identify an appropriate reference. 
        % In Networks: An Introduction, Newman (2010) cites:
        % Bearman et al. (2004) Chains of affection: the structure of
        % adolescent romantic and sexual networks
        % The webcite
        % https://web.archive.org/web/20210115130726/http://moreno.ss.uci.edu/data.html#adhealth
        % references:
        % Moody, James, "Peer influence groups: identifying dense clusters
        % in large networks," Social Networks, 2001, 23: 261-283. 
        % From 
        % 'The ADD HEALTH data are constructed from the in-school
        % questionnaire; 90,118 students representing 84 communities took
        % this survey in 1994-95. Some communities had only one school;
        % others had two. Where there are two schools in a community
        % students from one school were allowed to name friends in the
        % other, the "sister school."     
        % 
        % 'Each student was given a paper-and-pencil questionnaire and a
        % copy of a roster listing every student in the school and, if the
        % community had two schools, the student s provided with the roster
        % of the "sister" school. The name generator asked about five male
        % and five female friends separately. The question was, "List your
        % closest (male/female) friends. List your best (male/female)
        % friend first, then your next best friend, and so on. (girls/boys)
        % may include (boys/girls) who are friends and (boy/girl)
        % friends."'
        % 
        networkFile = 'udry_adol_health_large_conn_comp.mat';
        load([networkFolder, '/', networkFile], 'GCell')
        iiG = 67; G = GCell{iiG};
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency;
        nameStr = 'School 67';
        nameStr2 = 'High school friendship (2)';
        
    case 18
        
        % % Attractive (or at least interesting) result
        % 
        % N = 987
        % https://networkrepository.com/bn-mouse-kasthuri-graph-v4.php
        %
        % Amunts et al. (2013) BigBrain: An Ultrahigh-Resolution 3D Human
        % Brain Model 
        % From
        % https://networkrepository.com/bn-mouse-kasthuri-graph-v4.php
        % Edges represent fiber tracts that connect one vertex to another.
        %
        networkFile = 'bn-mouse-kasthuri_graph_v4.edges';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d', 'HeaderLines', 0); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency; N = size(A, 1); A(eye(N) == 1) = 0; G = graph(A);
        nameStr = 'Mouse visual cortex';
        nameStr2 = 'Visual cortex: Mouse';
        
    case 19
        
        % % Attractive result
        %
        % N = 1,458
        % https://networkrepository.com/bio-yeast-protein-inter.php
        %
        % Jeong et al. (2001), Lethality and centrality in protein
        % networks:
        % "The S. cerevisiae protein–protein interaction network we
        % investigate has 1,870 proteins as nodes, connected by 2,240
        % identified direct physical interactions, and is derived from
        % combined, non-overlapping data^3,4, obtained mostly by systematic
        % two-hybrid analyses^3."    
        %
        networkFile = 'bio-yeast-protein-inter.edges';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d', 'HeaderLines', 0); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency; N = size(A, 1); A(eye(N) == 1) = 0; G = graph(A);
        nameStr = 'Yeast protein';
        nameStr2 = 'Protein interaction: Yeast';
        
    case 20
        
        % % Attractive result
        %
        % N = 483
        % https://networkrepository.com/bio-DM-LC.php
        % 
        % Cho et al. (2013) WormNet v3: a network-assisted
        % hypothesis-generating server for Caenorhabditis elegans:
        % From 
        % http://www.inetbio.org/wormnet/downloadnetwork.php
        % "DM-LC	fly(D.melanogaster)	download	1129	Inferred Links by
        % small/medium-scale protein-protein interactions (collected from
        % protein-protein interaction data bases)"  
        % 
        networkFile = 'bio-DM-LC.edges';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d %f', 'HeaderLines', 0); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency; N = size(A, 1); A(eye(N) == 1) = 0; G = graph(A);
        nameStr = 'WormNet DM-LC';
        nameStr2 = 'Protein interaction: Fly (2)';
        
    case 21
        
        % % Attractive result
        %
        % N = 2,831
        % https://networkrepository.com/bio-DM-HT.php
        %
        % Cho et al. (2013) WormNet v3: a network-assisted
        % hypothesis-generating server for Caenorhabditis elegans:
        % From 
        % http://www.inetbio.org/wormnet/downloadnetwork.php
        % "DM-HT	fly(D.melanogaster)	download	4660	Inferred Links
        % by high-throughput protein-protein interactions"   
        %
        networkFile = 'bio-DM-HT.edges';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d %f', 'HeaderLines', 0); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency; N = size(A, 1); A(eye(N) == 1) = 0; G = graph(A);
        nameStr = 'WormNet DM-HT';
        nameStr2 = 'Protein interaction: Fly';
        
    case 22
        
        % % Attractive result
        %
        % N = 516
        % https://networkrepository.com/bio-diseasome.php
        % 
        % Goh et al. (2007) "In the "human disease network" (HDN) nodes represent
        % disorders, and two disorders are connected to each other if they
        % share at least one gene in which mutations are associated with
        % both disorders"  
        networkFile = 'bio-diseasome.mtx';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d', 'HeaderLines', 2); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency; N = size(A, 1); A(eye(N) == 1) = 0; G = graph(A);
        nameStr = 'Human diseases';
        nameStr2 = 'Human disease';
        
    case 23
        
        % % Attractive result
        %
        % N = 993
        % https://networkrepository.com/bio-CE-LC.php
        %
        % Cho et al. (2013) WormNet v3: a network-assisted
        % hypothesis-generating server for Caenorhabditis elegans:
        % From 
        % http://www.inetbio.org/wormnet/downloadnetwork.php
        % "CE-LC	worm(C.elegans)	download	1648	Inferred Links by
        % small/medium-scale protein-protein interactions (collected from
        % protein-protein interaction data bases)"
        %
        networkFile = 'bio-CE-LC.edges';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d', 'HeaderLines', 2); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency; N = size(A, 1); A(eye(N) == 1) = 0; G = graph(A);
        nameStr = 'WormNet CE-LC';
        nameStr2 = 'Protein interaction: Worm';
        
    case 24
        
        % % Attractive result
        %
        % N = 2,194
        % https://networkrepository.com/bio-CE-HT.php
        %
        % Cho et al. (2013) WormNet v3: a network-assisted
        % hypothesis-generating server for Caenorhabditis elegans:
        % From 
        % http://www.inetbio.org/wormnet/downloadnetwork.php
        % "CE-HT	worm(C.elegans)	download	2985	Inferred Links by
        % high-throughput protein-protein interactions"
        %
        networkFile = 'bio-CE-HT.edges';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d', 'HeaderLines', 2); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency; N = size(A, 1); A(eye(N) == 1) = 0; G = graph(A);
        nameStr = 'WormNet CE-HT';
        nameStr2 = 'Protein interaction: Worm (2)';
        
    case 25
        
        % % Attractive result
        % 
        % Pages linking to www.epa.gov
        % http://vlado.fmf.uni-lj.si/pub/networks/data/web/Epa.net
        % N = 4253
        %
        % http://vlado.fmf.uni-lj.si/pub/networks/data/web/Epa.net:
        % "This graph was constructed by expanding a 200-page response set
        % to a search engine query, as in the hub/authority algorithm. 
        % from Jon Kleinberg:
        %   http://www.cs.cornell.edu/courses/cs685/2002fa/
        % adapted for Pajek, V. Batagelj, March 19, 2006"  
        %
        networkFile = 'Epa.net';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, ' %d %d', 'HeaderLines', 4781); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency;
        nameStr = 'EPA';
        nameStr2 = 'Pages linking to www.epa.gov';
        
    case 26
        % % Attractive result
        %
        % Roads, Minnesota, Network Data Repository
        % https://networkrepository.com/road-minnesota.php
        % N = 2,640
        networkFile = 'road-minnesota.mtx';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d', 'HeaderLines', 15); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency;
        nameStr = 'Road network, Minnesota';
        nameStr2 = 'Road network: Minnesota';
        
    case 27
        % % Attractive result
        %
        % Roads, Europe, Network Data Repository
        % https://networkrepository.com/road-euroroad.php
        % N = 1,039
        networkFile = 'road-euroroad.edges';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d', 'HeaderLines', 2); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency;
        nameStr = 'Road network, Europe';
        nameStr2 = 'Road network: Europe';
        
    case 28
        % % Attractive result
        %
        % s208_st, Uri Alon
        % https://www.weizmann.ac.il/mcb/UriAlon/sites/mcb.UriAlon/files/uploads/CollectionsOfComplexNetwroks/s208_st.txt
        % N = 122
        networkFile = 's208_st.txt';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d %d', 'HeaderLines', 0); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency;
        nameStr = 'Circuit s208';
        nameStr2 = 'Circuit';
        
    case 29
        % % Attractive result
        %
        % s420_st, Uri Alon
        % https://www.weizmann.ac.il/mcb/UriAlon/sites/mcb.UriAlon/files/uploads/CollectionsOfComplexNetwroks/s420_st.txt
        % N = 252
        networkFile = 's420_st.txt';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d %d', 'HeaderLines', 0); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency;
        nameStr = 'Circuit s420';
        nameStr2 = 'Circuit (2)';
        
    case 30
        % % Attractive result
        %
        % s838_st, Uri Alon
        % https://www.weizmann.ac.il/mcb/UriAlon/sites/mcb.UriAlon/files/uploads/CollectionsOfComplexNetwroks/s838_st.txt
        % N = 512
        networkFile = 's838_st.txt';
        fileID = fopen([networkFolder, '/', networkFile]); C = textscan(fileID, '%d %d %d', 'HeaderLines', 0); fclose(fileID); from_list = C{:, 1}; to_list = C{:, 2}; minID = min([1, min(from_list), min(to_list)]);
        from_list = from_list + 1 - minID; to_list = to_list + 1 - minID; G = graph(from_list, to_list);
        [bin, binSize] = conncomp(G); iDsLCC = (binSize(bin) == max(binSize)); G = subgraph(G, iDsLCC);
        A = G.adjacency;
        nameStr = 'Circuit s838';
        nameStr2 = 'Circuit (3)';
end

end