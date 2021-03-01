import React from "react";
import { interpolateStructures } from "./utils/api";
import { LoadingPage } from "./components/LoadingPage";
import { StructureCardGrid } from "./components/StructureCardGrid";
import { StructureImage } from "./components/StructureImage";
import { ICollection, IInterpolatedParticle } from "./interfaces";
import { HorizontalCollapse } from "./components/HorizontalCollapse";
import { FormWrapper } from "./components/form";

export function InterpolationPage({
  structures,
  setStructures,
  collection,
  setCollection,
}: {
  structures?: string[];
  setStructures?(structures: string[]): void;
  collection?: ICollection<IInterpolatedParticle>;
  setCollection(collection: ICollection<IInterpolatedParticle>): void;
}) {
  const [optionsCollapsed, setOptionsCollapsed] = React.useState<boolean>(false);
  const [name, setName] = React.useState<string>("Interpolated");
  const [loading, setLoading] = React.useState<boolean>(false);

  return (
    <>
      <HorizontalCollapse
        label="Options"
        position="left"
        size="col-md-2"
        collapsed={optionsCollapsed}
        setCollapsed={setOptionsCollapsed}
      >
        <FormWrapper
          open={true}
          title="Use the continuous chemical space to interpolate between multiple structures."
          loading={loading}
          setLoading={setLoading}
          onSubmit={async () => {
            let data: IInterpolatedParticle[] = [];
            const validStructures = structures?.filter(Boolean);
            if (validStructures && validStructures.length > 0) {
              setLoading(true);
              data = await interpolateStructures(validStructures, 20);
            }
            setCollection({
              data,
              name,
            });
            setLoading(false);
          }}
        >
          {/* We rely on interpolation collection to be called "Interpolated" 
           <div className="form-group">
            <label htmlFor="interpolationNameInput">Name</label>
            <input
              type="text"
              className="form-control form-control-sm"
              id="interpolationNameInput"
              required
              value={name}
              onChange={(e) => setName(e.currentTarget.value)}
            />
          </div> */}
          <div className="form-group">
            <label htmlFor="inpolationStructures">Interpolate between</label>
            <textarea
              className="form-control form-control-sm"
              id="inpolationStructures"
              rows={3}
              onChange={(e) => setStructures?.(e.currentTarget.value.split("\n"))}
              value={structures?.join("\n") || ""}
            />
          </div>
          <div className="text-right">
            <button className="btn btn-primary" type="submit" disabled={loading}>
              {loading ? (
                <>
                  <span className="spinner-grow spinner-grow-sm" role="status" aria-hidden="true"></span> Loading...
                </>
              ) : (
                <>Go</>
              )}
            </button>
          </div>
          <div className="text-center">
            {structures && structures.length > 0
              ? structures
                  .map<React.ReactNode>((structure, i) => (
                    <StructureImage key={i} structure={structure} width="70px" height="70px" />
                  ))
                  .reduce((prev, curr) => [
                    prev,
                    <i className="fas  fa-fw fa-long-arrow-alt-right ml-2 mr-2"></i>,
                    curr,
                  ])
              : null}
          </div>
        </FormWrapper>
      </HorizontalCollapse>
      <div
        // className="col-md-10"
        style={{
          display: "flex",
          flex: 1,
          flexDirection: "column",
          overflow: "auto",
          marginLeft: 33,
        }}
      >
        <LoadingPage loading={loading} fallback="Please select structures for interpolation">
          {collection && collection.data.length > 0 ? (
            <StructureCardGrid
              collections={[collection]}
              tableClass="main-ranking"
              structureCardProps={(structure) => ({
                className: structure.scaffold ? "border-primary" : "",
              })}
            />
          ) : null}
        </LoadingPage>
      </div>
    </>
  );
}
