import * as React from "react";
import castArray from "lodash.castarray";
import Select from "react-select/creatable";

export interface GenericPlotSelectOption<T> {
  option: T | null | undefined;
  setOption(option: T | null): void;
}

export interface IPlotSelectExtension {
  component: React.FunctionComponent<{
    availableNearestNeighbors: string[];
    availableClusters: string[];
    availableProperties: string[];
  } & GenericPlotSelectOption<string>>;
  additionalOptions: {label: string, value: string}[];
}

export function PlotSelect<IsMulti extends boolean = false>({
  inline = true,
  multi,
  option,
  setOption,
  label,
  options,
  id,
  disabled,
  availableNearestNeighbors,
  availableClusters,
  availableProperties,
  extensions,
}: {
  inline?: boolean;
  options: string[];
  disabled?: boolean;
  id: string;
  label?: string;
  multi: IsMulti;
  availableNearestNeighbors: string[];
  availableClusters: string[];
  availableProperties: string[];
  extensions?: IPlotSelectExtension[];
} & (IsMulti extends true ? GenericPlotSelectOption<string[]> : GenericPlotSelectOption<string>)) {
  const values: string[] = React.useMemo(() => castArray(option).filter(Boolean) as string[], [option]);

  return (
    <div className={`d-flex position-relative ${inline ? 'col-sm-3' : 'col-12'}`}>
      {label ? (
        <label
          className={`my-1 me-2 align-self-center ${inline ? "" : "col-sm-2 d-flex"}`}
          style={{ whiteSpace: "nowrap" }}
          htmlFor={id}
        >
          {label}
        </label>
      ) : null}
      <div className={`my-1 me-sm-2 flex-fill ${inline ? "" : "col-sm-10"}`} style={{ minWidth: 100 }}>
        <div className="row">
          <div className="col">
            <Select<{ label: string; value: string }, IsMulti>
              menuPosition="fixed"
              isMulti={multi}
              isClearable
              styles={{ menu: (base) => ({ ...base, zIndex: 999 }) }}
              onCreateOption={(value) => {
                setOption(value as any);
              }}
              value={values.map((p) => ({
                label: p,
                value: p,
              }))}
              onChange={(e) => {
                if (!e) {
                  setOption(null);
                  return;
                }
                // @ts-ignore
                const values: { label: string; value: string }[] = castArray(e);
                if (multi) {
                  setOption(values.map((p) => p.value) as any);
                } else {
                  setOption(values?.[0]?.value as any);
                }
              }}
              options={[
                ...(extensions || []).map((e) => e.additionalOptions || []).flat(),
                ...options.map((option) => ({
                  label: option,
                  value: option,
                })),
              ]}
            />
          </div>
          {extensions?.map((e) => (
            <e.component
              setOption={(e) => {
                if(!e) {
                  return null;
                }
                const values: string[] = castArray(e);
                if (multi) {
                  setOption(values.map((p) => p) as any);
                } else {
                  setOption(values?.[0] as any);
                }
              }}
              availableNearestNeighbors={availableNearestNeighbors}
              availableClusters={availableClusters}
              availableProperties={availableProperties}
              option={values?.[0]}
            />
          ))}
        </div>
      </div>
    </div>
  );
}
