function model_rdm = get_model_rdm(categories_vec, cat_names, type)
% type:
%   'single-category' - within category dist=0, between 1
%   'low-level'       - within category=0
%                       between faces and watches=1/2
%                       between objects (non-watch) and animals=1/2
%                       between the rest=1
%   'semantic'        - within category=0
%                       between faces and animals=1/2
%                       between objects and watches=1/2
%                       between the rest=1
%   'face-vs-rest'    - within category=0
%                       between faces and anything else=1
%                       between non-face categories=1/2

if ~all(strcmp(cat_names,{'face','watch','object','animal'}))
    error('This function assumes your categories are: face,watch,object,animal (in this order)')
end
type = lower(type);
category_rdm = ones(4,4);
switch type
    case 'single-category' % nothing more to do
    case 'low-level'
        category_rdm([1,2],[1,2]) = 1/2;
        category_rdm([3,4],[3,4]) = 1/2;
    case 'semantic'
        category_rdm([1,4],[1,4]) = 1/2;
        category_rdm([2,3],[2,3]) = 1/2;
    case 'face-vs-rest'
        category_rdm(2:4,2:4) = 1/2;
end
category_rdm(~~eye(4,4)) = 0;
model_rdm = category_rdm(categories_vec, categories_vec);
end